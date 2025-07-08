from FeatureOutputWriter import FeatureFileWriter
from collections import deque, defaultdict
from dataclasses import dataclass, field
from typing import Optional, Generator
from Bio.SeqFeature import SeqFeature
from MultiCounter import MultiCounter
from Bio.SeqRecord import SeqRecord
from SiteFilter import SiteFilter
from utils import SiteVariantData
from BCBio import GFF
from site_variant_readers import (
    RNAVariantReader,
    Reditools2Reader,
    Reditools3Reader,
    Jacusa2Reader,
)
import progressbar
import argparse
import math

@dataclass(slots=True,frozen=True)
class QueueActionList:
    """
    Contains lists of features to "activate" (add the feature to an active_features dict) and features to "deactivate" (remove it from the active_features dict)
    """
    activate: list[SeqFeature] = field(default_factory=list)
    deactivate: list[SeqFeature] = field(default_factory=list)

class CountingContext():
    def __init__(
            self,
            feature_writer: FeatureFileWriter,
            filter: SiteFilter,
            use_progress_bar: bool
            ):
        self.active_features: dict[str,SeqFeature] = dict()
        self.feature_writer: FeatureFileWriter = feature_writer
        self.counters: defaultdict[str, MultiCounter] = defaultdict(self.default_counter_factory)
        self.use_progress_bar: bool = use_progress_bar
        self.action_queue: deque[tuple[int,QueueActionList]] = deque()
        self.filter = filter

        self.progbar: Optional[progressbar.ProgressBar] = None
        
        if self.use_progress_bar:
            self.progbar_increment = self._active_progbar_increment
        else:
            self.progbar_increment = self._inactive_progbar_increment

        return None
    
    def set_record(self, record: SeqRecord) -> None:
        assert len(self.active_features) == 0
        assert len(self.action_queue) == 0
        self.counters.clear()

        self.record: SeqRecord = record

        # Map positions to activation feature and deactivation actions
        position_actions: defaultdict[int,QueueActionList] = defaultdict(QueueActionList)
        for feature in record.features:
            self.load_action_queue(position_actions, feature, 1)

        # Create a queue of actions sorted by positions 
        self.action_queue.extend(sorted(position_actions.items()))

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
    
    def load_action_queue(self, location_actions: dict[int, QueueActionList], root_feature: SeqFeature, level: int) -> None:
        """
        Traverse a hierarchy stemming from a `root_feature`: Each visited feature is added to activation and deactivation actions in the `action_queue` according to the feature's `start` and `end` positions.
        """

        # Iterate the `parts` of a location for compatibility with `SimpleLocation` and `CompoundLocation`
        feature_strand: int = root_feature.location.parts[0].strand

        root_feature.level = level

        for part in root_feature.location.parts:
            if feature_strand != part.strand:
                raise Exception(f"feature {root_feature.id} contains parts on different strands ({feature_strand} and {part.strand}). I cannot work with this!")
            
            actions: QueueActionList = location_actions[int(part.start)]
            actions.activate.append(root_feature)

            actions: QueueActionList = location_actions[int(part.end)]
            actions.deactivate.append(root_feature)

        # Visit children
        for child in root_feature.sub_features:
            self.load_action_queue(location_actions, child, level + 1)

        return None

    def update_queues(self, new_position: int) -> None:
        """
        Activate or deactivate features depending on the actions in the `actions_queue` up to the current position.
        """

        visited_positions: int = 0
        while len(self.action_queue) > 0 and self.action_queue[0][0] <= new_position:
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
    
    def checkout(self, feature: SeqFeature) -> None:
        self.active_features.pop(feature.id, None)

        counter: Optional[MultiCounter] = self.counters.get(feature.id, None)

        if counter:
            self.feature_writer.write_feature_with_data(self.record, feature, counter)
            del self.counters[feature.id]
        else:
            self.feature_writer.write_feature_without_data(self.record, feature)

        for child in feature.sub_features:
            self.checkout(child)

        return None

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
        svdata: Optional[SiteVariantData] = reader.read()

        while svdata:
            self.update_queues(svdata.position)
            self.update_active_counters(svdata)
            svdata: SiteVariantData = reader.read()

        if not self.is_finished():
            last_position: int = self.action_queue[-1][0]
            self.update_queues(last_position)

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
    feature_output_filename: str = args.output + ".features.tsv" if args.output else "features.tsv"
    
    with (
        open(args.gff) as gff_handle,
        open(args.sites) as sv_handle,
        open(feature_output_filename, "w") as feature_output_handle
    ):
        match args.format:
            case "reditools2":
                sv_reader: RNAVariantReader = Reditools2Reader(sv_handle)
            case "reditools3":
                sv_reader: RNAVariantReader = Reditools3Reader(sv_handle)
            case "jacusa2":
                sv_reader: RNAVariantReader = Jacusa2Reader(sv_handle)
            case _:
                raise Exception(f'Unimplemented format "{args.format}"')
            
        feature_writer: FeatureFileWriter = FeatureFileWriter(feature_output_handle)
        feature_writer.write_comment(f"Input format: {args.format}")
        feature_writer.write_header()

        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        global_filter: SiteFilter = SiteFilter(
            cov_threshold=args.cov,
            edit_threshold=args.edit_threshold
        )

        ctx = CountingContext(feature_writer, global_filter, args.progress)
        for record in records:
            # manager: CountingContext = CountingContext(record, feature_writer, global_filter, args.progress)
            ctx.set_record(record)
            ctx.launch_counting(sv_reader)
