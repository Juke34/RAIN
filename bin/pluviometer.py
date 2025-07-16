#!/usr/bin/env python
from FeatureOutputWriter import FeatureFileWriter, AggregateFileWriter
from collections import deque, defaultdict
from dataclasses import dataclass, field
# from Bio.SeqFeature import SeqFeature
from typing import Optional, Generator
from MultiCounter import MultiCounter
from Bio.SeqRecord import SeqRecord
from SiteFilter import SiteFilter
from utils import SiteVariantData, location_union
from site_variant_readers import (
    RNAVariantReader,
    Reditools2Reader,
    Reditools3Reader,
    Jacusa2Reader,
)
from BCBio import GFF
from SeqFeature_extensions import SeqFeature
import progressbar
import argparse
import logging
import math

logger = logging.getLogger(__name__)

@dataclass(slots=True,frozen=True)
class QueueActionList:
    """
    Contains lists of features to "activate" (add the feature to an active_features dict) and features to "deactivate" (remove it from the active_features dict)
    """
    activate: list[SeqFeature] = field(default_factory=list)
    deactivate: list[SeqFeature] = field(default_factory=list)

def has_children_of_type(self: SeqFeature, target_type: str) -> None:
    """
    Return `True` if the feature contains sub-features of a specific target type
    """
    for child in self.sub_features:
        if child.type == target_type:
            return True
        
    return False



class CountingContext():
    def __init__(
            self,
            feature_writer: FeatureFileWriter,
            aggregate_writer: AggregateFileWriter,
            filter: SiteFilter,
            use_progress_bar: bool
            ):
        self.active_features: dict[str,SeqFeature] = dict()
        self.feature_writer: FeatureFileWriter = feature_writer
        self.aggregate_writer: AggregateFileWriter = aggregate_writer
        self.counters: defaultdict[str, MultiCounter] = defaultdict(self.default_counter_factory)
        self.aggregate_counters: defaultdict[str, MultiCounter] = defaultdict(self.default_counter_factory)
        self.use_progress_bar: bool = use_progress_bar
        self.action_queue: deque[tuple[int,QueueActionList]] = deque()
        self.filter: SiteFilter = filter
        self.svdata: Optional[SiteVariantData] = None
        self.deactivation_list: list[SeqFeature] = []

        self.progbar: Optional[progressbar.ProgressBar] = None
        
        if self.use_progress_bar:
            self.progbar_increment = self._active_progbar_increment
        else:
            self.progbar_increment = self._inactive_progbar_increment

        return None
    
    def set_record(self, record: SeqRecord) -> None:
        logging.info(f"Switching to record {record.id}")
        if len(self.active_features) != 0:
            logging.info(f"Features in active queque not cleared prior to switching!\n{','.join(self.active_features.keys())}")
            self.active_features.clear()
        if len(self.action_queue) != 0:
            logging.info(f"Positions in action queque not cleared prior to switching!\n{','.join(str(x[0]) for x in self.action_queue)}")
            self.action_queue.clear()
        self.counters.clear()

        self.record: SeqRecord = record

        # Map positions to activation feature and deactivation actions
        logging.info("Pre-processing features")
        position_actions: defaultdict[int,QueueActionList] = defaultdict(QueueActionList)
        for feature in record.features:
            self.load_action_queue(position_actions, feature, 1, ['.'])

        # Create a queue of actions sorted by positions 
        self.action_queue.extend(sorted(position_actions.items()))
        logging.info(f"Pre-processing complete. Action positions detected: {len(self.action_queue)}")

        if self.use_progress_bar:
            max_value: int = len(self.action_queue)
            self.progbar = progressbar.ProgressBar(
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
    
    def default_counter_factory(self) -> MultiCounter:
        return MultiCounter(self.filter)
    
    def _active_progbar_increment(self, i: int) -> None:
        """Just a wrapper for the `ProgressBar.increment` method"""
        self.progbar.increment(i)

        return None
    
    def _inactive_progbar_increment(self, i:int) -> None:
        """Dummy method that does nothing when the progress bar is deactivated"""
        pass
    
    def load_action_queue(self, location_actions: dict[int, QueueActionList], root_feature: SeqFeature, level: int, parent_list: list[str]) -> None:
        """
        Traverse a hierarchy stemming from a `root_feature`: Each visited feature is added to activation and deactivation actions in the `action_queue` according to the feature's `start` and `end` positions.
        """

        # Iterate over the `parts` of a location for compatibility with `SimpleLocation` and `CompoundLocation`
        feature_strand: int = root_feature.location.parts[0].strand

        root_feature.level = level
        root_feature.parent_list = parent_list

        # if "chimaera" not in root_feature.type:
        if level == 1:
            logging.info(f"Creating chimaera of feature {root_feature.id}")
            root_feature.make_chimaera()

        for part in root_feature.location.parts:
            if feature_strand != part.strand:
                raise Exception(f"feature {root_feature.id} contains parts on different strands ({feature_strand} and {part.strand}). I cannot work with this!")
            
            actions: QueueActionList = location_actions[int(part.start)]
            actions.activate.append(root_feature)

            actions: QueueActionList = location_actions[int(part.end)]
            actions.deactivate.append(root_feature)

        # Visit children
        if hasattr(root_feature, "sub_features"):
            for child in root_feature.sub_features:
                self.load_action_queue(location_actions, child, level + 1, parent_list + [root_feature.id])

        return None
    
    def state_update_cycle(self, new_position: int) -> None:
        """
        Activate or deactivate features depending on the actions in the `actions_queue` up to the current position.
        """

        visited_positions: int = 0

        while len(self.action_queue) > 0 and self.action_queue[0][0] < new_position:    # Use < instead of <= because of Python's right-exclusive indexing
            _, actions = self.action_queue.popleft()
            visited_positions += 1
            
            for feature in actions.activate:
                self.active_features[feature.id] = feature

            for feature in actions.deactivate:
                if feature.level == 1:
                    self.checkout(feature)

                self.active_features.pop(feature.id, None)

        self.progbar_increment(visited_positions)

        return None
    
    def flush_queues(self) -> None:
        # last_position: int = len(self.record)
        logging.info(f"Actions remaining in action queue: {len(self.action_queue)}. Flushing action queue")
        while len(self.action_queue) > 0:
            _, actions = self.action_queue.popleft()

            for feature in actions.activate:
                self.active_features[feature.id] = feature

            for feature in actions.deactivate:
                if feature.level == 1:
                    self.checkout(feature)

                self.active_features.pop(feature.id, None)
        logging.info(f"Actions remaining in action queue: {len(self.action_queue)} after flushing action queue")

        return None
    
    def checkout(self, feature: SeqFeature) -> defaultdict[str,MultiCounter]:
        self.active_features.pop(feature.id, None)

        # Counter for the feature itself
        counter: Optional[MultiCounter] = self.counters.get(feature.id, None)

        if counter:
            self.feature_writer.write_row_with_data(self.record.id, feature, counter)
            del self.counters[feature.id]
        else:
            self.feature_writer.write_row_without_data(self.record.id, feature)

        # Aggregation counters from the feature's sub-features
        if feature.level == 1:
            aggregation_counters: dict[str,MultiCounter] = self.aggregate_level1(feature)
        else:
            aggregation_counters: dict[str,MultiCounter] = self.aggregate_children(feature)

        # Recursively check-out children
        for child in feature.sub_features:
            self.checkout(child)

        return aggregation_counters
    
    def aggregate_level1(self, feature: SeqFeature) -> dict[str,MultiCounter]:
        aggregation_counters: defaultdict[str,MultiCounter] = defaultdict(self.default_counter_factory)

        # List of tuples of transcript-like sub-features. In each tuple:
        # - 0: ID of the sub-feature
        # - 1: Type of the sub-feature
        # - 2: Total length of the sub-feature
        transcript_like_children:list[tuple[str,str,int]] = feature.get_transcript_like()   # Custom method added to the class

        # Select the transcript-like feature that is representative of this gene.
        # If there are CDS sub-features, select the onte with greatest total CDS length. Elsewise, select the sub-feature with the greatest total exon length.
        representative_feature_id: str = ""
        has_cds: bool = False
        max_total_length: bool = 0

        for child_id, child_type, child_length in transcript_like_children:
            if child_type == "cds":
                if has_cds:
                    if child_length > max_total_length:
                        representative_feature_id = child_id
                        max_total_length = child_length
                else:
                    representative_feature_id = child_id
                    max_total_length = child_length
                    has_cds = True
            elif child_type == "exon":
                if has_cds:
                    continue
                elif child_length > max_total_length:
                        representative_feature_id = child_id
                        max_total_length = child_length

        # Perform aggregations, selecting only the "representative feature"
        for child in feature.sub_features:
            # Compute aggregates in the child. Recursively aggregates on all its children.
            aggregation_counters_from_child = self.aggregate_children(child)

            if child.id == representative_feature_id:
                # Merge the aggregates from the child with all the other aggregates under this feature
                for child_aggregation_type, child_aggregation_counter in aggregation_counters_from_child.items():
                    aggregation_counter: MultiCounter = aggregation_counters[child_aggregation_type]
                    aggregation_counter.merge(child_aggregation_counter)

        #     self.aggregate_writer.write_row_with_data(
        #         self.record.id,
        #         feature,
        #         aggregation_type,
        #         "",
        #         aggregation_counter
        #     )

        self.aggregate_writer.write_row_with_data(self.record.id, feature, aggregation_counters)

        return aggregation_counters


    def aggregate_children(self, feature: SeqFeature) -> dict[str,MultiCounter]:
        aggregation_counters: defaultdict[str,MultiCounter] = defaultdict(self.default_counter_factory)

        for child in feature.sub_features:
            # Compute aggregates in the child
            aggregation_counters_from_child = self.aggregate_children(child)

            # Merge the aggregates from the child with all the other aggregates under this feature
            for child_aggregation_type, child_aggregation_counter in aggregation_counters_from_child.items():
                aggregation_counter: MultiCounter = aggregation_counters[child_aggregation_type]
                aggregation_counter.merge(child_aggregation_counter)

            # Add the child itself to the aggregations
            aggregation_counter: MultiCounter = aggregation_counters[child.type]
            feature_counter: Optional[MultiCounter] = self.counters.get(child.id, None)
            if feature_counter:
                aggregation_counter.merge(feature_counter)
            
        self.aggregate_writer.write_row_with_data(self.record.id, feature, aggregation_counters)
        # for aggregation_type, aggregation_counter in aggregation_counters.items():
        #     self.aggregate_writer.write_row_with_data(
        #         self.record.id,
        #         feature,
        #         aggregation_type,
        #         "",
        #         aggregation_counter
        #     )

        return aggregation_counters


    def update_active_counters(self, site_data: SiteVariantData) -> None:
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
    
    def launch_counting(self, reader: RNAVariantReader) -> None:
        # self.svdata: Optional[SiteVariantData] = self.svdata if self.svdata else reader.read() 
        next_svdata: Optional[SiteVariantData] = reader.seek_record(self.record.id)
        if next_svdata:
            logging.info(f"Found a site variant data matching the record {self.record.id}")

        while next_svdata and next_svdata.seqid == self.record.id:
            self.svdata = next_svdata
            self.state_update_cycle(self.svdata.position)
            self.update_active_counters(self.svdata)
            next_svdata: Optional[SiteVariantData] = reader.read()

        if self.svdata:
            logging.info(f"Last position in record {self.record.id}: {self.svdata.position}")

        if not self.is_finished():
            # last_position: int = len(self.record)
            # logging.info(f"Updating queues up to position {last_position}")
            # last_position: int = max(max(map(lambda x: x.location.end, self.active_features.values())), max(map(lambda x: x[0], self.action_queue[0])))
            # self.update_queues(last_position)
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
        "-t",
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
        "--progress", action="store_true", default=False, help="Display progress bar"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args: argparse.Namespace = parse_cli_input()
    LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    # logging.basicConfig(format=LOGGING_FORMAT)
    log_filename: str = args.output + ".pluviometer.log" if args.output else "pluviometer.log"
    logging.basicConfig(filename=log_filename, level=logging.INFO, format=LOGGING_FORMAT)
    logging.info(f"Pluviometer started. Log file: {log_filename}")
    feature_output_filename: str = args.output + ".features.tsv" if args.output else "features.tsv"
    aggregate_output_filename: str = args.output + ".aggregates.tsv" if args.output else "aggregates.tsv"
    
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
                # sv_reader: RNAVariantReader = Reditools3Reader(sv_handle)
            case "jacusa2":
                reader_factory = Jacusa2Reader
                # sv_reader: RNAVariantReader = Jacusa2Reader(sv_handle)
            case _:
                raise Exception(f'Unimplemented format "{args.format}"')
            
        feature_writer: FeatureFileWriter = FeatureFileWriter(feature_output_handle)
        feature_writer.write_comment(f"Input format: {args.format}")
        feature_writer.write_header()

        aggregate_writer = AggregateFileWriter(aggregate_output_handle)
        aggregate_writer.write_comment(f"Input format: {args.format}")
        aggregate_writer.write_header()

        logging.info("Load GFF3 file")
        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        global_filter: SiteFilter = SiteFilter(
            cov_threshold=args.cov,
            edit_threshold=args.edit_threshold
        )

        ctx = CountingContext(feature_writer, aggregate_writer, global_filter, args.progress)
        logging.info("Fetching records...")
        for i, record in enumerate(records):
            with open(args.sites) as sv_handle:
                sv_reader = reader_factory(sv_handle)
                logging.info(f"Record {record.id} setup")
                ctx.set_record(record)
                logging.info(f"Start counting on record {record.id}")
                ctx.launch_counting(sv_reader)
                logging.info(f"Ended counting on record {record.id}")
                
