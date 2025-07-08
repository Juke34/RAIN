from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from typing import TextIO, Optional, Generator
import numpy as np
from numpy.typing import NDArray
import argparse
from MultiCounter import MultiCounter
from SiteFilter import SiteFilter
from utils import SiteVariantData, condense
from contextlib import nullcontext
import sys
from site_variant_readers import (
    RNAVariantReader,
    Reditools2Reader,
    Reditools3Reader,
    Jacusa2Reader,
)
from FeatureOutputWriter import RainFileWriter
from collections import deque, defaultdict
import progressbar
import math

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
        help="Name of the output file (leave empty to write to stdout)",
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
        "--progress", action="store_true", default="false", help="Display progress bar"
    )

    return parser.parse_args()


class RecordManager:
    def __init__(self, record: SeqRecord, writer: RainFileWriter, filter: SiteFilter):

        self.record: SeqRecord = record
        self.writer: RainFileWriter = writer
        """Dict of multicounters with feature ids as keys"""
        self.counters: defaultdict[str, MultiCounter] = defaultdict(self.counter_factory)
        """Set of features whose multicounters are currently being updated"""
        self.active_features: dict[str,SeqFeature] = dict()
        self.filter: SiteFilter = filter

        # Flatten the feature hierarchy
        feature_list: list[SeqFeature] = []

        for feature in record.features:
            self._flatten_hierarchy(feature_list, feature, 1)

        self.nb_targets = len(feature_list)

        nb_targets_d_format = math.floor(math.log(self.nb_targets, 10))

        # Create deques of features sorted by their start position and their end position
        # These deques are used for loading and unloading features from the `active_features` set
        feature_list.sort(key=lambda feature: feature.location.start)
        self.activation_deque: deque[tuple[int, list[SeqFeature]]] = condense(feature_list, "start")
        print(f"Features to activate: {len(self.activation_deque)}")

        feature_list.sort(key=lambda feature: feature.location.end)
        self.deactivation_deque: deque[tuple[int, list[SeqFeature]]] = condense(feature_list, "end")
        print(f"Features to deactivate: {len(self.deactivation_deque)}")


        # self.features_start_first: deque[SeqFeature] = deque(
        #     sorted(feature_list, key=lambda feature: feature.location.start),
        # )

        # self.features_end_first: deque[SeqFeature] = deque(
        #     sorted(feature_list, key=lambda feature: feature.location.end),
        # )

        # self.switchlist: deque[tuple[int,str,str]] = deque()


        self.use_progress_bar = args.progress
        if self.use_progress_bar:
            self.progress_bar: progressbar.ProgressBar = progressbar.ProgressBar(
                max_value=self.nb_targets,
                widgets=[
                    f"Record {self.record.id}: Processed target feature ",
                    progressbar.Counter(
                        format=f"%(value)0{nb_targets_d_format}d out of %(max_value)d"
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

    def counter_factory(self) -> MultiCounter:
        return MultiCounter(self.filter)

    def _flatten_hierarchy(self, feature_list, feature, level) -> None:
        """Add elements to a flattened list of features by preorder traversal of a feature hierarchy"""
        feature_list.append(feature)
        feature.level = level
        for child in feature.sub_features:
            self._flatten_hierarchy(feature_list, child, level + 1)

        return None
    
    def update_active_counters(self, site_data: SiteVariantData) -> None:
        """
        Update the multicounters matching the ID of features in the `active_features` set.
        A new multicounter is created if no matching ID is found.
        """
        for feature_key, feature in self.active_features.items():
            counter = self.counters[feature_key]
            # if counter:
            #     print(f"counter found for feature {key}")
            if feature.location.strand == site_data.strand:
                counter.update(site_data)

        return None
    
    # def update_queues(self, new_position: int) -> None:
    #     remove_list: list[str] = []

    #     delete_list: list[str] = []
    #     for key, feature in self.active_features.items():
    #         if feature.location.end < new_position:
    #             delete_list.append(key)
        
    #     for key in delete_list:
    #         del self.active_features[key]
    #         # self.progress_bar.next()

    #     while len(self.features_start_first) > 0 and self.features_start_first[0].location.start <= new_position:
    #         new_feature: SeqFeature = self.features_start_first.popleft()

    #         if new_feature.location.end > new_position:
    #             self.active_features[new_feature.id] = new_feature
    #         else:
    #             # self.progress_bar.next()
    #             pass
    #             # print(f"activated feature {new_feature.id}")

    #     # while len(self.features_end_first) > 0 and self.features_end_first[0].location.end < new_position:
    #     #     old_feature: SeqFeature = self.features_end_first.popleft()

    #     #     remove_list.append(old_feature.id)
    #     #     # self.active_features.discard(old_feature.id)
    #     #     if old_feature.level == 1:
    #     #         # print(f"deactivated feature {old_feature.id}")
    #     #         self.checkout(old_feature)

    #     # for elem in remove_list:
    #     #     self.active_features.discard(elem)

    #     return None

    def update_queues(self, new_position: int) -> None:
        while len(self.activation_deque) > 0 and self.activation_deque[0][0] <= new_position:
            feature_list = self.activation_deque.popleft()[1]
            for feature in feature_list:
                if feature.location.end < new_position:
                    if args.progress:
                        self.progress_bar.next()
                else:
                    self.active_features[feature.id] = feature

        while len(self.deactivation_deque) > 0 and self.deactivation_deque[0][0] < new_position:
            feature_list = self.deactivation_deque.popleft()[1]
            for feature in feature_list:
                self.active_features.pop(feature.id, None)
                if feature.level == 1:
                    self.checkout(feature)
                # if args.progress:
                #     self.progress_bar.next()


    def checkout(self, feature: SeqFeature) -> None:
        # print(f"checking out feature {feature.id}")
        if feature.id in self.active_features:
            del self.active_features[feature.id]
            if args.progress:
                self.progress_bar.next()

        counter: Optional[MultiCounter] = self.counters.get(feature.id, MultiCounter(self.filter))
        # print(self.counters.keys())
        if counter:
            # print(f"writing data for feature {feature.id}")
            self.writer.write_feature_data(self.record, feature, counter)
            del self.counters[feature.id]

        for child in feature.sub_features:
            self.checkout(child)

        return None
    
    def is_finished(self) -> bool:
        result = len(self.active_features) + len(self.activation_deque) + len(self.deactivation_deque) == 0
        # if len(self.activation_deque) == 0:
        #     print("Activation queue empty")
        # if len(self.deactivation_deque) == 0:
        #     print("Deactivation queue empty")
        # if result:
        #     print("All queues are empty")
        
        return result

    def launch_counting(self, reader: RNAVariantReader) -> None:
        svdata: Optional[SiteVariantData] = reader.read()

        while svdata and not self.is_finished():
            # print(self.active_features)
            self.update_queues(svdata.position)
            self.update_active_counters(svdata)
            svdata: SiteVariantData = reader.read()
            # print(f"Active counters: {len(self.counters)} - Active features: {len(self.active_features)}")



if __name__ == "__main__":
    args = parse_cli_input()
    with (
        open(args.gff) as gff_handle,
        open(args.sites) as sv_handle,
        open(args.output, "w")
        if len(args.output) > 0
        else nullcontext(sys.stdout) as output_handle,
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

        writer: RainFileWriter = RainFileWriter(output_handle)
        writer.write_comment(f"input format: {args.format}")
        writer.write_header()

        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        global_filter: SiteFilter = SiteFilter(
            cov_threshold=args.cov, edit_threshold=args.edit_threshold
        )

        for record in records:
            manager: RecordManager = RecordManager(record, writer, global_filter)
            manager.launch_counting(sv_reader)

        print(f"Active: {len(manager.active_features)}")
        # if len(manager.active_features) > 0:
            # print(*map(lambda x: x.id, manager.active_features))
        print(f"To activate: {len(manager.activation_deque)}")
        # if len(manager.activation_deque) > 0:
            # print(*map(lambda x: x.id, manager.activation_deque))
        print(f"To deactivate: {len(manager.deactivation_deque)}")
        # if len(manager.deactivation_deque) > 0:
            # print(*map(lambda x: x.id, manager.deactivation_deque))
