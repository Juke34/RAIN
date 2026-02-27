#!/usr/bin/env python
from .rain_file_writers import FeatureFileWriter, AggregateFileWriter
from typing import Any, Optional, Generator, Callable, TextIO
from Bio.SeqFeature import SimpleLocation, CompoundLocation
from .SeqFeature_extensions import SeqFeature
from collections import deque, defaultdict
from dataclasses import dataclass, field
from .multi_counter import MultiCounter
from Bio.SeqRecord import SeqRecord
from .site_filter import SiteFilter
from .rna_site_variant_readers import (
    RNASiteVariantReader,
    Reditools2Reader,
    Reditools3Reader,
    Jacusa2Reader,
)
from .utils import RNASiteVariantData
from natsort import natsorted
import multiprocessing
from BCBio import GFF
from os import remove
import progressbar
import tempfile
import argparse
import logging
import math

logger = logging.getLogger(__name__)


@dataclass(slots=True, frozen=True)
class QueueActionList:
    """
    Contains lists of features to "activate" (add the feature to an active_features dict) and features to "deactivate" (remove it from the active_features dict)
    """

    activate: list[SeqFeature] = field(default_factory=list)
    deactivate: list[SeqFeature] = field(default_factory=list)


@dataclass(slots=True)
class AggregatePositions:
    """
    Contains the extreme positions (start, end, strand) for an aggregate of features
    """
    start: int = float('inf')  # Will be replaced by minimum start position
    end: int = 0  # Will be replaced by maximum end position
    strand: int = 0  # Strand from parent feature
    
    def update_from_feature(self, feature: SeqFeature) -> None:
        """Update positions with a feature's location"""
        if hasattr(feature.location, 'start') and hasattr(feature.location, 'end'):
            self.start = min(self.start, int(feature.location.start))
            self.end = max(self.end, int(feature.location.end))
            if self.strand == 0 and hasattr(feature.location, 'strand'):
                self.strand = feature.location.strand if feature.location.strand else 0
    
    def merge(self, other: 'AggregatePositions') -> None:
        """Merge with another AggregatePositions"""
        if other.start != float('inf'):
            self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)
        if self.strand == 0:
            self.strand = other.strand
    
    def to_strings(self) -> tuple[str, str, str]:
        """Convert to string representation for output"""
        start_str = str(self.start + 1) if self.start != float('inf') else "."  # Convert to 1-based
        end_str = str(self.end) if self.end > 0 else "."
        strand_str = str(self.strand) if self.strand != 0 else "."
        return start_str, end_str, strand_str


def merge_aggregation_counter_dicts(
    dst: defaultdict[str, MultiCounter], src: defaultdict[str, MultiCounter]
) -> None:
    for src_type, src_counter in src.items():
        dst_counter: MultiCounter = dst[src_type]
        dst_counter.merge(src_counter)

    return None


class RecordCountingContext:
    """
    Handles counting at record (~ chromosome) level.
    """

    def __init__(
        self,
        feature_writer: FeatureFileWriter,
        aggregate_writer: AggregateFileWriter,
        filter: SiteFilter,
        use_progress_bar: bool,
    ):
        self.aggregate_writer: AggregateFileWriter = aggregate_writer
        """Aggregate counting output is written to a temporary file through this object"""

        self.feature_writer: FeatureFileWriter = feature_writer
        """Feature counting output is written to a temporary file through this object"""

        self.filter: SiteFilter = filter
        """Filter applied to the reads. Filters contain a mutable state, so they should not be shared between concurrent RecordCountingContext instances"""

        self.longest_isoform_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for the record-level longest isoform aggregates. Keys correspond to aggregate feature types"""

        self.chimaera_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for the record-level chimaeric aggregates. Keys correspond to aggregate feature types"""

        self.all_isoforms_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for the record-level all isoforms aggregates. Keys correspond to aggregate feature types"""

        # New: Aggregate counters by (parent_type, aggregate_type) for intermediate-level aggregation
        self.longest_isoform_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for longest isoform aggregates grouped by parent type. Keys are (parent_type, aggregate_type) tuples"""

        self.chimaera_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for chimaeric aggregates grouped by parent type. Keys are (parent_type, aggregate_type) tuples"""

        self.all_isoforms_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """Dictionary of counters for all isoforms aggregates grouped by parent type. Keys are (parent_type, aggregate_type) tuples"""

        self.total_counter: MultiCounter = MultiCounter(self.filter)
        """The counter that keeps tallies of editions on all the sites of the record, regardless of feature"""

        self.active_features: dict[str, SeqFeature] = dict()
        """This dictonary holds all the features whose counters are must be updated at the current genomic position. Keys are feature IDs"""

        self.counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        """This dictionary contains the feature-level counters. New counters are created on the fly as needed by the use of a default dictionary. Keys are feature IDs and entries can be deleted after checkout to save memory"""

        self.use_progress_bar: bool = use_progress_bar
        """Toggle for progress bar. This is more useful during development"""

        self.action_queue: deque[tuple[int, QueueActionList]] = deque()
        """This queue (implemented with a deque) holds a list of all the activating or deactivating actions to execute at the genomic positions downstream of the current position"""

        self.svdata: Optional[RNASiteVariantData] = None
        """The site variant data object holding the information about variants mapped to the current genomic position"""


        self.progbar_increment: Callable = (
            self._active_progbar_increment
            if self.use_progress_bar
            else self._inactive_progbar_increment
        )
        """The method to be used to update the progress bar. When no progress bar is used, a dummy method is called"""

        return None

    def update_aggregate_counters(
        self, aggregate_counter_tag: str, new_counters: defaultdict[str, MultiCounter]
    ) -> None:
        """
        Update a dict of aggregate counters by merging them with matching items in another dict of counters.
        The tag refers to the different kinds of aggregate counters: "all_isoforms", "chimaera", and "longest_isoform"

        Arguments:
            aggregate_counter_tag: str  Prefix for selecting one of the aggregate counter dictionaries of the different aggregation modes: "all_isoforms", "longest_isoform", and "chimaera"
            new_counters: defaultdict[str, MultiCounter]    Defaultdict of counters to merge with the counters selected by the aggregate counter tag.
        """
        aggregate_counters: defaultdict[str, MultiCounter] = getattr(
            self, aggregate_counter_tag + "_aggregate_counters"
        )
        for counter_type, new_counter in new_counters.items():
            target_counter: MultiCounter = aggregate_counters[counter_type]
            target_counter.merge(new_counter)

        return None

    def set_record(self, record: SeqRecord) -> None:
        """
        Assign a specific record to the RecordCountingContext.
        If the context was already assigned to another record, the queues and counters are cleared before proceeding (this was used in a previous version that recycled the record counting context, but it isn't used since parallelization was implemented).
        Makes the reader seek the first entry matching the record ID, and initiates the pre-processing of the features.
        1. Parse all the features in the record and load the resulting actions into the action queue
        2. 
        """
        logging.info(f"Switching to record {record.id}")
        if len(self.active_features) != 0:
            logging.info(
                f"Features in active queque not cleared prior to switching!\n{','.join(self.active_features.keys())}"
            )
            self.active_features.clear()
        if len(self.action_queue) != 0:
            logging.info(
                f"Positions in action queque not cleared prior to switching!\n{','.join(str(x[0]) for x in self.action_queue)}"
            )
            self.action_queue.clear()
        self.counters.clear()

        self.record: SeqRecord = record

        # Map positions to activation feature and deactivation actions
        logging.info("Pre-processing features")
        location_actions: defaultdict[int, QueueActionList] = defaultdict(QueueActionList)
        for feature in record.features:
            # Loading is recursive
            self.load_location_actions_dict(
                location_actions=location_actions,
                root_feature=feature, 
                level=1,
                parent_list=["."]   # "." represents the record root, i.e. the parent of a level-1 position
                )

        # Create a queue of actions sorted by positions
        self.action_queue.extend(sorted(location_actions.items()))
        logging.info(
            f"Pre-processing complete. Action positions detected: {len(self.action_queue)}"
        )

        if self.use_progress_bar:
            max_value: int = len(self.action_queue)
            self.progbar: progressbar.ProgressBar = progressbar.ProgressBar(
                max_value=max_value,
                widgets=[
                    f"Record {self.record.id}: Position ",
                    progressbar.Counter(
                        format=f"%(value)0{math.floor(math.log(max_value, 10))}d out of %(max_value)d"
                    ),
                    " (",
                    progressbar.Percentage(),
                    ") ",
                    progressbar.Bar("█", "|", "|"),
                    " ",
                    progressbar.Timer(),
                    " - ",
                    progressbar.SmoothingETA(),
                ],
                poll_interval=1,  # Updates every 1 second
            )

        return None

    def _active_progbar_increment(self, i: int) -> None:
        """Just a wrapper for the `ProgressBar.increment` method"""
        self.progbar: progressbar.ProgressBar
        self.progbar.increment(i)

        return None

    def _inactive_progbar_increment(self, i: int) -> None:
        """Dummy method that does nothing when the progress bar is deactivated"""
        pass

    def load_location_actions_dict(
        self,
        location_actions: dict[int, QueueActionList],
        root_feature: SeqFeature,
        level: int,
        parent_list: list[str],
    ) -> None:
        """
        Traverse a hierarchy stemming from a `root_feature`: Each visited feature is added to activation and deactivation actions in the `action_queue` according to the feature's `start` and `end` positions.
        """

        # Iterate over the `parts` of a location for compatibility with `SimpleLocation` and `CompoundLocation`
        assert root_feature.location
        feature_strand: Optional[int] = root_feature.location.parts[0].strand

        root_feature.level = level
        if not hasattr(root_feature, "is_chimaera"):
            root_feature.is_chimaera = False
        root_feature.parent_list = parent_list

        if level == 1:
            root_feature.make_chimaeras(self.record.id)

        for part in root_feature.location.parts:
            if feature_strand != part.strand:
                raise Exception(
                    f"feature {root_feature.id} contains parts on different strands ({feature_strand} and {part.strand}). I cannot work with this!"
                )

        old_part: Optional[SimpleLocation | CompoundLocation] = None

        for part in root_feature.location.parts:
            if old_part:
                if old_part.contains(part.start) or old_part.contains(part.stop):
                    raise Exception(
                        f"feature {root_feature.id} has a compound location containing overlapping parts. There must be no overlapping."
                    )

            actions: QueueActionList = location_actions[int(part.start)]
            actions.activate.append(root_feature)

            actions = location_actions[int(part.end)]
            actions.deactivate.append(root_feature)

        # Visit children
        if hasattr(root_feature, "sub_features"):
            for child in root_feature.sub_features:
                self.load_location_actions_dict(
                    location_actions, child, level + 1, parent_list + [root_feature.id]
                )

        return None

    def state_update_cycle(self, new_position: int) -> None:
        """
        Activate or deactivate features depending on the actions in the `actions_queue` up to the current position.
        """

        visited_positions: int = 0

        while (
            len(self.action_queue) > 0 and self.action_queue[0][0] < new_position
        ):  # Use < instead of <= because of Python's right-exclusive indfgexing
            _, actions = self.action_queue.popleft()
            visited_positions += 1

            for feature in actions.activate:
                self.active_features[feature.id] = feature

            for feature in actions.deactivate:
                if feature.level == 1:
                    self.checkout(feature, None)

                self.active_features.pop(feature.id, None)

        self.progbar_increment(visited_positions)

        return None

    def flush_queues(self) -> None:
        """
        Finish processing all the queues without checking for more data from the RNA variant input file
        """

        logging.info(
            f"Actions remaining in action queue: {len(self.action_queue)}. Flushing action queue"
        )
        while len(self.action_queue) > 0:
            _, actions = self.action_queue.popleft()

            for feature in actions.activate:
                self.active_features[feature.id] = feature

            for feature in actions.deactivate:
                if feature.level == 1:
                    self.checkout(feature, None)

                self.active_features.pop(feature.id, None)
        logging.info(
            f"Actions remaining in action queue: {len(self.action_queue)} after flushing action queue"
        )

        return None

    def checkout(self, feature: SeqFeature, parent_feature: Optional[SeqFeature]) -> None:
        """
        Deactivate and post-process a feature (and its sub-features, recursively) and write to output the associated counting data.
        
        Checkout should be triggered after a level-1 feature is completely cleared. Aggregation of the level-1 feature and all its descentants proceeds.
        """

        # Deactivate the feature
        self.active_features.pop(feature.id, None)

        # Counter for the feature itself
        feature_counter: Optional[MultiCounter] = self.counters.get(feature.id, None)

        assert self.record.id

        if feature_counter:
            if feature.is_chimaera:
                assert parent_feature  # A chimaera must always have a parent feature (a gene)
                self.aggregate_writer.write_row_chimaera_with_data(
                    self.record.id, feature, parent_feature, feature_counter
                )
                self.chimaera_aggregate_counters[feature.type].merge(feature_counter)
                # Also track by parent_type
                self.chimaera_aggregate_counters_by_parent_type[(parent_feature.type, feature.type)].merge(feature_counter)
            else:
                self.feature_writer.write_row_with_data(self.record.id, feature, feature_counter)
            del self.counters[feature.id]
        else:
            if feature.is_chimaera:
                assert parent_feature
                self.aggregate_writer.write_row_chimaera_without_data(
                    self.record.id, feature, parent_feature
                )
            else:
                self.feature_writer.write_row_without_data(self.record.id, feature)

        # all_isoforms_aggregation_counters: Optional[defaultdict[str, MultiCounter]] = None

        assert self.record.id  # Placate Pylance

        # Aggregation counters from the feature's sub-features
        if feature.level == 1:
            (
                level1_longest_isoform_aggregation_counters,
                level1_all_isoforms_aggregation_counters,
                level1_longest_isoform_aggregation_positions,
                level1_all_isoforms_aggregation_positions,
            ) = self.aggregate_level1(feature)
            merge_aggregation_counter_dicts(
                self.all_isoforms_aggregate_counters, level1_all_isoforms_aggregation_counters
            )

            # Also merge into by-parent-type dictionaries
            parent_type = feature.type
            for aggregate_type, aggregate_counter in level1_longest_isoform_aggregation_counters.items():
                self.longest_isoform_aggregate_counters_by_parent_type[(parent_type, aggregate_type)].merge(aggregate_counter)
            for aggregate_type, aggregate_counter in level1_all_isoforms_aggregation_counters.items():
                self.all_isoforms_aggregate_counters_by_parent_type[(parent_type, aggregate_type)].merge(aggregate_counter)

            for (
                aggregate_type,
                aggregate_counter,
            ) in level1_longest_isoform_aggregation_counters.items():
                # Get positions for this aggregate type
                start_str, end_str, strand_str = ".", ".", "."
                if aggregate_type in level1_longest_isoform_aggregation_positions:
                    start_str, end_str, strand_str = level1_longest_isoform_aggregation_positions[aggregate_type].to_strings()
                
                self.aggregate_writer.write_metadata(
                    seq_id=self.record.id,
                    parent_ids=".," + feature.id,
                    aggregate_id=feature.longest_isoform + "-longest_isoform",
                    parent_type=feature.type,
                    aggregate_type=aggregate_type,
                    aggregation_mode="longest_isoform",
                    start=start_str,
                    end=end_str,
                    strand=strand_str,
                )
                self.aggregate_writer.write_counter_data(aggregate_counter)

            for (
                aggregate_type,
                aggregate_counter,
            ) in level1_all_isoforms_aggregation_counters.items():
                # Get positions for this aggregate type
                start_str, end_str, strand_str = ".", ".", "."
                if aggregate_type in level1_all_isoforms_aggregation_positions:
                    start_str, end_str, strand_str = level1_all_isoforms_aggregation_positions[aggregate_type].to_strings()
                
                self.aggregate_writer.write_metadata(
                    seq_id=self.record.id,
                    parent_ids=".," + feature.id,
                    aggregate_id=feature.id + "-all_isoforms",
                    parent_type=feature.type,
                    aggregate_type=aggregate_type,
                    aggregation_mode="all_isoforms",
                    start=start_str,
                    end=end_str,
                    strand=strand_str,
                )
                self.aggregate_writer.write_counter_data(aggregate_counter)
        else:
            feature_aggregation_counters, feature_aggregation_positions = self.aggregate_children(feature)
            for aggregate_type, aggregate_counter in feature_aggregation_counters.items():
                # Get positions for this aggregate type
                start_str, end_str, strand_str = ".", ".", "."
                if aggregate_type in feature_aggregation_positions:
                    start_str, end_str, strand_str = feature_aggregation_positions[aggregate_type].to_strings()
                
                self.aggregate_writer.write_metadata(
                    seq_id=self.record.id,
                    parent_ids=",".join(feature.parent_list),
                    aggregate_id=feature.id,
                    parent_type=parent_feature.type,
                    aggregate_type=aggregate_type,
                    aggregation_mode="feature",
                    start=start_str,
                    end=end_str,
                    strand=strand_str,
                )
                self.aggregate_writer.write_counter_data(aggregate_counter)

        # Recursively check-out children
        for child in feature.sub_features:
            self.checkout(child, feature)

        return None

    def aggregate_level1(
        self, feature: SeqFeature
    ) -> tuple[defaultdict[str, MultiCounter], defaultdict[str, MultiCounter], dict[str, AggregatePositions], dict[str, AggregatePositions]]:
        """
        Perform aggregation for level-1 feature.
        
        It handles the special cases caused by the different modes of aggregation (all_isoforms, longest_isoform, and chimaera).
        Triggers the recursive aggregation of the feature descendants (features of higher levels)
        Returns: (longest_isoform_counters, all_isoforms_counters, longest_isoform_positions, all_isoforms_positions)
        """

        level1_longest_isoform_aggregation_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        level1_all_isoforms_aggregation_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        level1_longest_isoform_aggregation_positions: dict[str, AggregatePositions] = {}
        level1_all_isoforms_aggregation_positions: dict[str, AggregatePositions] = {}

        # List of tuples of transcript-like sub-features. In each tuple:
        # - 0: ID of the sub-feature
        # - 1: Type of the sub-feature
        # - 2: Total length of the sub-feature
        transcript_like_children: list[tuple[str, str, int]] = (
            feature.get_transcript_like()
        )  # Custom method added to the class

        # Select the transcript-like feature that is representative of this gene: the longest isoform.
        # If there are CDS sub-features, select the onte with greatest total CDS length. Elsewise, select the sub-feature with the greatest total exon length.
        longest_isoform_id: str = ""
        has_cds: bool = False
        max_total_length: int = 0

        for child_id, child_type, child_length in transcript_like_children:
            if child_type == "CDS":
                if has_cds:
                    if child_length > max_total_length:
                        longest_isoform_id = child_id
                        max_total_length = child_length
                else:
                    longest_isoform_id = child_id
                    max_total_length = child_length
                    has_cds = True
            elif child_type == "exon":
                if has_cds:
                    continue
                elif child_length > max_total_length:
                    longest_isoform_id = child_id
                    max_total_length = child_length

        feature.longest_isoform = longest_isoform_id

        logging.info(
            f"Record {self.record.id}, gene {feature.id}: The longest isoform is {longest_isoform_id}."
        )

        # Perform aggregations
        for child in feature.sub_features:
            # Compute aggregates in the child. Recursively aggregates on all its children.
            aggregation_counters_from_child, aggregation_positions_from_child = (
                self.aggregate_children(child)
            )

            merge_aggregation_counter_dicts(
                level1_all_isoforms_aggregation_counters, aggregation_counters_from_child
            )
            
            # Merge positions for all_isoforms
            for aggregate_type, positions in aggregation_positions_from_child.items():
                if aggregate_type not in level1_all_isoforms_aggregation_positions:
                    level1_all_isoforms_aggregation_positions[aggregate_type] = AggregatePositions(strand=feature.location.strand if hasattr(feature.location, 'strand') else 0)
                level1_all_isoforms_aggregation_positions[aggregate_type].merge(positions)

            if child.id == longest_isoform_id:
                # Merge the aggregates from the child with all the other aggregates under this feature
                merge_aggregation_counter_dicts(
                    level1_longest_isoform_aggregation_counters, aggregation_counters_from_child
                )
                # Merge positions for longest_isoform
                for aggregate_type, positions in aggregation_positions_from_child.items():
                    if aggregate_type not in level1_longest_isoform_aggregation_positions:
                        level1_longest_isoform_aggregation_positions[aggregate_type] = AggregatePositions(strand=feature.location.strand if hasattr(feature.location, 'strand') else 0)
                    level1_longest_isoform_aggregation_positions[aggregate_type].merge(positions)

        # Merge the feature-level aggregation counters into the record-level aggregation counters
        self.update_aggregate_counters(
            "longest_isoform", level1_longest_isoform_aggregation_counters
        )

        return (
            level1_longest_isoform_aggregation_counters,
            level1_all_isoforms_aggregation_counters,
            level1_longest_isoform_aggregation_positions,
            level1_all_isoforms_aggregation_positions,
        )

    def aggregate_children(self, feature: SeqFeature) -> tuple[defaultdict[str, MultiCounter], dict[str, AggregatePositions]]:
        """
        Perform aggregation of features. Not intended to handle level-1 features (see `aggregate_level1`)
        Returns a tuple of (counters_dict, positions_dict) where positions_dict contains extreme positions for each aggregate type
        """

        aggregation_counters: defaultdict[str, MultiCounter] = defaultdict(
            DefaultMultiCounterFactory(self.filter)
        )
        aggregation_positions: dict[str, AggregatePositions] = {}

        for child in feature.sub_features:
            # Compute aggregates in the child
            aggregation_counters_from_child, aggregation_positions_from_child = self.aggregate_children(child)

            # Merge the aggregates from the child with all the other aggregates under this feature
            for (
                child_aggregation_type,
                child_aggregation_counter,
            ) in aggregation_counters_from_child.items():
                aggregation_counter: MultiCounter = aggregation_counters[child_aggregation_type]
                aggregation_counter.merge(child_aggregation_counter)
                
                # Merge positions
                if child_aggregation_type not in aggregation_positions:
                    aggregation_positions[child_aggregation_type] = AggregatePositions(strand=feature.location.strand if hasattr(feature.location, 'strand') else 0)
                if child_aggregation_type in aggregation_positions_from_child:
                    aggregation_positions[child_aggregation_type].merge(aggregation_positions_from_child[child_aggregation_type])

            # Add the child itself to the aggregations
            aggregation_counter: MultiCounter = aggregation_counters[child.type]
            feature_counter: Optional[MultiCounter] = self.counters.get(child.id, None)
            if feature_counter:
                aggregation_counter.merge(feature_counter)
                
            # Update positions with the child feature's position
            if child.type not in aggregation_positions:
                aggregation_positions[child.type] = AggregatePositions(strand=feature.location.strand if hasattr(feature.location, 'strand') else 0)
            aggregation_positions[child.type].update_from_feature(child)

        return aggregation_counters, aggregation_positions

    def update_active_counters(self, site_data: RNASiteVariantData) -> None:
        """
        Update the multicounters matching the ID of features in the `active_features` set.
        A new multicounter is created if no matching ID is found.
        """
        for feature_key, feature in self.active_features.items():
            if feature.location.strand == site_data.strand:
                counter: MultiCounter = self.counters[feature_key]
                counter.update(site_data)

        return None

    def is_finished(self) -> bool:
        return len(self.action_queue) > 0 or len(self.active_features) > 0

    def launch_counting(self, reader: RNASiteVariantReader) -> None:
        next_svdata: Optional[RNASiteVariantData] = reader.seek_record(self.record.id)
        if next_svdata:
            logging.info(f"Record {self.record.id} · Found site variant data matching the record")

        while next_svdata and next_svdata.seqid == self.record.id:
            self.svdata = next_svdata
            self.state_update_cycle(self.svdata.position)
            self.update_active_counters(self.svdata)
            self.total_counter.update(self.svdata)
            next_svdata: Optional[RNASiteVariantData] = reader.read()

        if self.svdata:
            logging.info(
                f"Record {self.record.id} · Last position in genome reached: {self.svdata.position}"
            )

        if not self.is_finished():
            self.flush_queues()

        if self.use_progress_bar:
            self.progbar.finish()

        return None


def parse_cli_input() -> argparse.Namespace:
    """Parse command line input"""

    parser = argparse.ArgumentParser(description="Rain counter")
    parser.add_argument(
        "--sites",
        "-s",
        type=str,
        required=True,
        help="File containing per-site base alteration data",
    )
    parser.add_argument(
        "--gff",
        "-g",
        type=str,
        required=True,
        help="Reference genome annotations (GFF3 file)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="",
        type=str,
        help="Prefix for the names of the output files",
    )
    parser.add_argument(
        "--format",
        "-f",
        type=str,
        choices=["reditools2", "reditools3", "jacusa2", "sapin"],
        default="reditools3",
        help="Sites file format",
    )
    parser.add_argument(
        "--cov",
        "-c",
        type=int,
        default=0,
        help="Site coverage threshold for counting editions",
    )
    parser.add_argument(
        "--edit_threshold",
        "-T",
        type=int,
        default=1,
        help="Minimum number of edited reads for counting a site as edited",
    )
    parser.add_argument(
        "--aggregation_mode",
        "-a",
        type=str,
        default="all",
        choices=["all", "cds_longest"],
        help='Mode for aggregating counts: "all" aggregates features of every transcript; "cds_longest" aggregates features of the longest CDS or non-coding transcript',
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads (actually, processes) to use for parallel computing",
    )
    parser.add_argument(
        "--progress", action="store_true", default=False, help="Display progress bar"
    )

    return parser.parse_args()


class DefaultMultiCounterFactory:
    """Callable class to enable the pickling of a MultiCounter factory for multiprocessing (lambdas cannot be pickled)"""

    def __init__(self, filter: SiteFilter):
        self.filter: SiteFilter = filter

        return None

    def __call__(self):
        return MultiCounter(self.filter)


# Global variables for multiprocessing
args = None
reader_factory = None


def init_worker(args_value, reader_factory_value):
    """Initialize global variables in each worker process"""
    global args, reader_factory
    args = args_value
    reader_factory = reader_factory_value


def run_job(record: SeqRecord) -> dict[str, Any]:
    """
    A wrapper function for performing counting parallelized by record. The return value is a dict containing all the information needed for integrating
    the output of all records after the computations are finished.
    """
    assert record.id  # Stupid assertion for pylance
    logging.info(f"Record {record.id} · Record parsed. Counting beings.")

    tmp_feature_output_file: str = tempfile.mkstemp()[1]
    tmp_aggregate_output_file: str = tempfile.mkstemp()[1]

    with (
        open(args.sites) as sv_handle,
        open(tmp_feature_output_file, "w") as tmp_feature_output_handle,
        open(tmp_aggregate_output_file, "w") as tmp_aggregate_output_handle,
    ):
        # Set up output
        feature_writer: FeatureFileWriter = FeatureFileWriter(tmp_feature_output_handle)
        aggregate_writer: AggregateFileWriter = AggregateFileWriter(tmp_aggregate_output_handle)

        # Set up context
        filter: SiteFilter = SiteFilter(cov_threshold=args.cov, edit_threshold=args.edit_threshold)

        record_ctx: RecordCountingContext = RecordCountingContext(
            feature_writer, aggregate_writer, filter, args.progress
        )

        # Count
        reader: RNASiteVariantReader = reader_factory(sv_handle)
        record_ctx.set_record(record)
        record_ctx.launch_counting(reader)

        # Write aggregate counter data of the record
        aggregate_writer.write_rows_with_data(
            record.id,
            ["."],
            ".",
            ".",
            "longest_isoform",
            record_ctx.longest_isoform_aggregate_counters,
        )
        aggregate_writer.write_rows_with_data(
            record.id,
            ["."],
            ".",
            ".",
            "all_isoforms",
            record_ctx.all_isoforms_aggregate_counters,
        )
        aggregate_writer.write_rows_with_data(
            record.id,
            ["."],
            ".",
            ".",
            "chimaera",
            record_ctx.chimaera_aggregate_counters,
        )

        # Write aggregate counter data by parent type
        aggregate_writer.write_rows_with_data_by_parent_type(
            record.id,
            ["."],
            ".",
            "longest_isoform",
            record_ctx.longest_isoform_aggregate_counters_by_parent_type,
        )
        aggregate_writer.write_rows_with_data_by_parent_type(
            record.id,
            ["."],
            ".",
            "all_isoforms",
            record_ctx.all_isoforms_aggregate_counters_by_parent_type,
        )
        aggregate_writer.write_rows_with_data_by_parent_type(
            record.id,
            ["."],
            ".",
            "chimaera",
            record_ctx.chimaera_aggregate_counters_by_parent_type,
        )

        # Write the total counter data of the record. A dummy dict needs to be created to use the `write_rows_with_data` method
        total_counter_dict: defaultdict[str, MultiCounter] = defaultdict(
            lambda: MultiCounter(filter)
        )
        total_counter_dict["."] = record_ctx.total_counter
        aggregate_writer.write_rows_with_data(
            record.id, ["."], ".", ".", "all_sites", total_counter_dict
        )

    return {
        "record_id": record.id,
        "tmp_feature_output_file": tmp_feature_output_file,
        "tmp_aggregate_output_file": tmp_aggregate_output_file,
        "chimaera_aggregate_counters": record_ctx.chimaera_aggregate_counters,
        "longest_isoform_aggregate_counters": record_ctx.longest_isoform_aggregate_counters,
        "all_isoforms_aggregate_counters": record_ctx.all_isoforms_aggregate_counters,
        "chimaera_aggregate_counters_by_parent_type": record_ctx.chimaera_aggregate_counters_by_parent_type,
        "longest_isoform_aggregate_counters_by_parent_type": record_ctx.longest_isoform_aggregate_counters_by_parent_type,
        "all_isoforms_aggregate_counters_by_parent_type": record_ctx.all_isoforms_aggregate_counters_by_parent_type,
        "total_counter": record_ctx.total_counter,
    }


def main():
    global args
    global reader_factory

    args = parse_cli_input()
    LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    log_filename: str = args.output + "_pluviometer.log" if args.output else "pluviometer.log"
    logging.basicConfig(filename=log_filename, level=logging.INFO, format=LOGGING_FORMAT)
    logging.info(f"Pluviometer started. Log file: {log_filename}")
    feature_output_filename: str = args.output + "_features.tsv" if args.output else "features.tsv"
    aggregate_output_filename: str = (
        args.output + "_aggregates.tsv" if args.output else "aggregates.tsv"
    )

    with (
        open(args.gff) as gff_handle,
        open(feature_output_filename, "w") as feature_output_handle,
        open(aggregate_output_filename, "w") as aggregate_output_handle,
    ):
        match args.format:
            case "reditools2":
                reader_factory = Reditools2Reader
            case "reditools3":
                reader_factory = Reditools3Reader
            case "jacusa2":
                reader_factory = Jacusa2Reader
            case _:
                raise Exception(f'Unimplemented format "{args.format}"')

        logging.info("Parsing GFF3 file...")
        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        genome_filter: SiteFilter = SiteFilter(
            cov_threshold=args.cov, edit_threshold=args.edit_threshold
        )
        genome_total_counter: MultiCounter = MultiCounter(genome_filter)

        genome_longest_isoform_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )
        genome_all_isoforms_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )
        genome_chimaera_aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )

        # By-parent-type genome-level aggregates
        genome_longest_isoform_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )
        genome_all_isoforms_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )
        genome_chimaera_aggregate_counters_by_parent_type: defaultdict[tuple[str, str], MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )

        feature_writer: FeatureFileWriter = FeatureFileWriter(feature_output_handle)
        aggregate_writer: AggregateFileWriter = AggregateFileWriter(aggregate_output_handle)

        feature_writer.write_header()
        aggregate_writer.write_header()

        with multiprocessing.Pool(processes=args.threads, initializer=init_worker, initargs=(args, reader_factory)) as pool:
            # Run the jobs and save the result of each record in a dict stored in the `record_data` list
            record_data_list: list[dict[str, Any]] = pool.map(run_job, records)

        # Sort record results in "natural" order
        record_data_list = natsorted(record_data_list, key=lambda x: x["record_id"])

    with (
        open(feature_output_filename, "a") as feature_output_handle,
        open(aggregate_output_filename, "a") as aggregate_output_handle,
    ):
        for record_data in record_data_list:
            logging.info(f"Record {record_data['record_id']} · Merging temporary output files...")

            with open(record_data["tmp_feature_output_file"]) as tmp_output_handle:
                feature_output_handle.write(tmp_output_handle.read())
            remove(record_data["tmp_feature_output_file"])

            with open(record_data["tmp_aggregate_output_file"]) as tmp_output_handle:
                aggregate_output_handle.write(tmp_output_handle.read())
            remove(record_data["tmp_aggregate_output_file"])

            # Update the genome's aggregate counters from the record data aggregate counters
            for record_aggregate_type, record_aggregate_counter in record_data[
                "longest_isoform_aggregate_counters"
            ].items():
                genome_aggregate_counter: MultiCounter = genome_longest_isoform_aggregate_counters[
                    record_aggregate_type
                ]
                genome_aggregate_counter.merge(record_aggregate_counter)

            for record_aggregate_type, record_aggregate_counter in record_data[
                "chimaera_aggregate_counters"
            ].items():
                genome_aggregate_counter: MultiCounter = genome_chimaera_aggregate_counters[
                    record_aggregate_type
                ]
                genome_aggregate_counter.merge(record_aggregate_counter)

            merge_aggregation_counter_dicts(
                genome_all_isoforms_aggregate_counters,
                record_data["all_isoforms_aggregate_counters"],
            )

            # Update the genome's by-parent-type aggregate counters
            for (parent_type, aggregate_type), record_aggregate_counter in record_data[
                "longest_isoform_aggregate_counters_by_parent_type"
            ].items():
                genome_longest_isoform_aggregate_counters_by_parent_type[(parent_type, aggregate_type)].merge(
                    record_aggregate_counter
                )

            for (parent_type, aggregate_type), record_aggregate_counter in record_data[
                "all_isoforms_aggregate_counters_by_parent_type"
            ].items():
                genome_all_isoforms_aggregate_counters_by_parent_type[(parent_type, aggregate_type)].merge(
                    record_aggregate_counter
                )

            for (parent_type, aggregate_type), record_aggregate_counter in record_data[
                "chimaera_aggregate_counters_by_parent_type"
            ].items():
                genome_chimaera_aggregate_counters_by_parent_type[(parent_type, aggregate_type)].merge(
                    record_aggregate_counter
                )

            # Update the genome's total counter from the record data total counter
            genome_total_counter.merge(record_data["total_counter"])

        logging.info("Writing genome totals...")

        aggregate_writer: AggregateFileWriter = AggregateFileWriter(aggregate_output_handle)

        # Write genomic counts
        aggregate_writer.write_rows_with_data(
            ".", ["."], ".", ".", "longest_isoform", genome_longest_isoform_aggregate_counters
        )
        aggregate_writer.write_rows_with_data(
            ".", ["."], ".", ".", "all_isoforms", genome_all_isoforms_aggregate_counters
        )
        aggregate_writer.write_rows_with_data(
            ".", ["."], ".", ".", "chimaera", genome_chimaera_aggregate_counters
        )

        # Write genomic counts by parent type
        aggregate_writer.write_rows_with_data_by_parent_type(
            ".", ["."], ".", "longest_isoform", genome_longest_isoform_aggregate_counters_by_parent_type
        )
        aggregate_writer.write_rows_with_data_by_parent_type(
            ".", ["."], ".", "all_isoforms", genome_all_isoforms_aggregate_counters_by_parent_type
        )
        aggregate_writer.write_rows_with_data_by_parent_type(
            ".", ["."], ".", "chimaera", genome_chimaera_aggregate_counters_by_parent_type
        )

        # Write the genomic total. A dummy dict needs to be created to use the `write_rows_with_data` method
        genomic_total_counter_dict: defaultdict[str, MultiCounter] = defaultdict(
            lambda: MultiCounter(genome_filter)
        )
        genomic_total_counter_dict["."] = genome_total_counter
        aggregate_writer.write_rows_with_data(
            ".", ["."], ".", ".", "all_sites", genomic_total_counter_dict
        )

        logging.info("Program finished")
if __name__ == "__main__":
        print ()"Program finished")
    main()