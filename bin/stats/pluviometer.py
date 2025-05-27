#!/usr/bin/env python

#####################
# Conventions:
#   - Both internally and in the output, when base type information is represented in an array, the array indices follow the order ACGT
#   - Always use explicit return
#   - Test that collections are empty (or not) by comparing their lengths to zero (don't test "if x" as though x where False or None)
#####################

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import numpy as np
from numpy.typing import NDArray
from typing import TextIO, NamedTuple, Optional, Generator
from site_variant_readers import RNAVariantReader, ReditoolsReader
import argparse
from contextlib import nullcontext
import sys
from utils import SiteVariantData, BASE_TYPES, EDIT_TYPES
from collections import deque
import progressbar
import math
from feature_aggregator import FeatureAggregator


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
        choices=["reditools", "jacusa2", "sapin"],
        default="reditools",
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
        help="Mode for aggregating counts: \"all\" aggregates features of every transcript; \"cds_longest\" aggregates features of the longest CDS or non-coding transcript",
    )

    return parser.parse_args()


class SiteFilter:
    def __init__(self, cov_threshold: int, edit_threshold: int) -> None:
        self.cov_threshold: int = cov_threshold
        self.edit_threshold: int = edit_threshold
        self.frequencies: NDArray[np.int32] = np.zeros(4, np.int32)

    def apply(self, variant_data: SiteVariantData) -> None:
        if variant_data.coverage >= self.cov_threshold:
            np.copyto(self.frequencies, variant_data.frequencies * variant_data.frequencies >= self.edit_threshold)
        else:
            self.frequencies.fill(0)

        return None


class MultiCounter:
    """Holds the counter data and logic for a feature, feature aggregate, or record"""

    def __init__(self, site_filter: SiteFilter) -> None:
        """
        Tallies of the number of times that a read has been mapped a site in the genome with base X.
        It is *not* the number of times a base appears in a read
        """
        self.ref_coverage_by_base_type: NDArray[np.int32] = np.zeros(4, dtype=np.int32)

        """
        Number of times that each base type occurs in the reference genome
        """
        self.ref_base_freqs: NDArray[np.int32] = np.zeros(4, dtype=np.int32)

        """
        Tallies of the numbers of reads per edit type
        This is a numpy matrix where the rows represent the reference base and the columns the edited base
        Rows and column indices correspond to bases in alphabetic order (ACGT)
        Row-columns corresponding to the same base (e.g. (0,0) -> (A,A)) do not represent edits, and should remain 0
        """
        self.edit_read_freqs: NDArray[np.int32] = np.zeros((4, 4), dtype=np.int32)

        self.edit_site_freqs: NDArray[np.int32] = np.zeros((4, 4), dtype=np.int32)

        self.filter = site_filter

        return None

    def update(self, variant_data: SiteVariantData) -> None:
        """Increment the counters from the data in a SiteVariantData object."""
        i: int = variant_data.reference

        # TODO: Decide whether to count the coverage of bases when it is below the coverage threshold used for counting edits
        self.ref_coverage_by_base_type[i] += variant_data.coverage
        self.ref_base_freqs[i] += variant_data.coverage > 0

        self.edit_read_freqs[i, :] += variant_data.frequencies

        self.filter.apply(variant_data)
        self.edit_site_freqs[i, :] += self.filter.frequencies

        return None

    def report(self, output_handle: TextIO) -> int:
        b = 0

        # Write the number of covered sites
        b += output_handle.write(str(self.ref_base_freqs.sum()))
        output_handle.write("\t")

        # Write reference base frequencies
        b += write_base_array(output_handle, self.ref_base_freqs)
        output_handle.write("\t")

        # Write edited sites
        b += write_edit_array(output_handle, self.edit_site_freqs)
        output_handle.write("\t")

        # Write sums of base coverages
        b += write_base_array(output_handle, self.ref_coverage_by_base_type)
        output_handle.write("\t")

        # Write edit frequencies
        b += write_edit_array(output_handle, self.edit_read_freqs)
        # output_handle.write("\t")

        # # Write proportion of edited reads
        # b += write_edit_array(output_handle, self.compute_proportions())

        return b

    def write_ref_coverage_by_base(self, output_handle: TextIO) -> int:
        """Writes out the reference base read frequencies in comma-separated format"""
        return output_handle.write(
            ",".join(str(value) for value in self.ref_coverage_by_base_type)
        )

    def write_edit_freqs(self, output_handle: TextIO) -> int:
        """Writes the edit frequencies in a comma-separated format"""
        # Yes, this looks like lisp
        return output_handle.write(
            ",".join(
                ",".join(
                    # Skip indices where i == j, because they don't represent editions
                    str(self.edit_read_freqs[i, j])
                    for j in filter(lambda j: j != i, range(4))
                )
                for i in range(4)
            )
        )

    def write_edited_sites(self, output_handle: TextIO) -> int:
        return output_handle.write(
            ",".join(str(value) for value in self.edit_site_freqs)
        )

    def compute_proportions(self) -> NDArray[np.float64]:
        return np.divide(
            self.edit_read_freqs,
            # Add ones to bases with zero reads to avoid division by 0
            self.ref_coverage_by_base_type + (self.ref_coverage_by_base_type == 0),
        )


def write_base_array(output_handle: TextIO, x: NDArray) -> int:
    return output_handle.write(",".join(str(value) for value in x))


def write_edit_array(output_handle: TextIO, x: NDArray) -> int:
    return output_handle.write(
        ",".join(
            ",".join(
                # Skip indices where i == j, because they don't represent editions
                str(x[i, j])
                for j in filter(lambda j: j != i, range(4))
            )
            for i in range(4)
        )
    )


class FeatureGroupManager:
    """
    Manages a collection of level-1 features and all their sub-features, plus the associated counters
    A feature collection contains roots
    """

    def __init__(self, root: SeqFeature, recman: "RecordManager"):
        """Level-1 features"""
        self.roots: list[SeqFeature] = [root]
        self.counters: dict[str, MultiCounter] = dict()
        self.recman: "RecordManager" = recman
        self.nb_remaining_target_nodes: int = 0

        self.recman.aggregator.aggregate_sub_features(root)

        self.init_children(root)

    def init_children(self, feature: SeqRecord) -> None:
        for child in feature.sub_features:
            self.init_children(child)

        if feature.location.strand:
            self.recman.sorted_target_features.append(ManagedFeature(feature, self))
            if feature.id in self.counters:
                counter = self.counters[feature.id]   # Doesn't do anything!
            else:
                counter = MultiCounter(self.recman.filter)
                self.counters[feature.id] = counter

            self.nb_remaining_target_nodes += 1

        return None

    def add_root(self, feature: SeqRecord) -> None:
        if self.roots[0].id != feature.id:
            raise Exception(
                f"Cannot add feature with id {feature.id} to the roots of feature group manager with id {self.roots[0].id}"
            )

        self.roots.append(feature)

        self.recman.aggregator.aggregate_sub_features(feature)

        self.init_children(feature)

        return None

    def update_counters(
        self, feature: SeqFeature, variant_data: SiteVariantData
    ) -> None:
        if feature.location.strand == variant_data.strand:
            self.counters[feature.id].update(variant_data)

        return None

    def write_feature_data(self, feature: SeqFeature, grandparent:str, parent: str, output_handle: TextIO) -> None:
        b: int = 0
        b += output_handle.write(
            f"{self.recman.record.id}\t{grandparent}\t{parent}\t{feature.id}\t{feature.type}\t"
            +
            # Shift start location to GFF 1-based index
            f"{feature.location.start + 1}\t{feature.location.end}\t{feature.location.strand}\t"
        )

        b += self.counters[feature.id].report(output_handle)
        b += output_handle.write("\n")

        return b

    def checkout(self, feature: SeqFeature, grandparent: str, parent: str, output_handle: TextIO) -> None:
        if feature.id in self.counters:
            self.write_feature_data(feature, grandparent, parent, output_handle)
            self.counters.pop(feature.id, None)

        for child in feature.sub_features:
            self.checkout(child, parent, feature.id, output_handle)

        return None

    def update_progress(self, output_handle: TextIO) -> None:
        self.nb_remaining_target_nodes -= 1
        if self.nb_remaining_target_nodes == 0:
            for root in self.roots:
                self.checkout(root, ".", ".", output_handle)

        return None


class ManagedFeature(NamedTuple):
    """Wrapper for keeping the association between a feature and its feature group manager"""

    feature: SeqFeature
    manager: FeatureGroupManager


class QueueUpdates(NamedTuple):
    activated: int
    dectivated: int
    skipped: int


class RecordManager:
    def __init__(self, record: SeqRecord, global_filter: SiteFilter, output_handle: TextIO, aggregation_mode: str):
        self.record: SeqRecord = record
        self.pos: int = 0
        self.final_pos: int = len(record.seq)
        self.downstream_queue: deque[ManagedFeature] = deque()
        self.active_queue: deque[ManagedFeature] = deque()
        self.next_start: int = -1
        self.next_end: int = -1
        self.sorted_target_features: list[ManagedFeature] = []
        self.feature_managers: dict[str, FeatureGroupManager] = dict()
        self.output_handle: TextIO = output_handle
        self.filter: SiteFilter = global_filter
        self.counter: MultiCounter = MultiCounter(self.filter)
        self.aggregator: FeatureAggregator = FeatureAggregator(aggregation_mode)

        self.start_pos = record.features[0].location.start

        for feature in record.features:
            self.start_pos = min(self.start_pos, feature.location.start)
            if feature.id in self.feature_managers:
                self.feature_managers[feature.id].add_root(feature)
            else:
                feature_manager = FeatureGroupManager(feature, self)
                self.feature_managers[feature.id] = feature_manager

        self.sorted_target_features.sort(
            key=lambda x: (x.feature.location.start, x.feature.location.end)
        )
        self.nb_targets = len(self.sorted_target_features)

        nb_targets_d_format = math.floor(math.log(self.nb_targets, 10))

        self.progress_bar = progressbar.ProgressBar(
            max_value=self.nb_targets,
            widgets=[
                "Processed target feature ",
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
            poll_interval=0.001,
        )

        self.downstream_queue.extend(self.sorted_target_features)

        self.next_start = self.downstream_queue[0].feature.location.start
        self.next_end = self.downstream_queue[0].feature.location.end

        return None

    def is_finished(self) -> bool:
        return len(self.downstream_queue) + len(self.active_queue) == 0

    def load_active_queue(self, new: ManagedFeature) -> None:
        for i, old in enumerate(self.active_queue):
            if new.feature.location.end < old.feature.location.end:
                self.active_queue.insert(i, new)

                return None

        self.active_queue.append(new)

        return None

    def update_queues(self, pos: int) -> QueueUpdates:
        activated: int = 0
        deactivated: int = 0
        skipped: int = 0

        while (
            len(self.active_queue) > 0
            and self.active_queue[0].feature.location.end < pos
        ):
            item: ManagedFeature = self.active_queue.popleft()
            deactivated += 1
            item.manager.update_progress(self.output_handle)
            self.progress_bar.next()

        while (
            len(self.downstream_queue) > 0
            and self.downstream_queue[0].feature.location.end < pos
        ):
            item = self.downstream_queue.popleft()
            skipped += 1
            item.manager.update_progress(self.output_handle)
            self.progress_bar.next()

        while (
            len(self.downstream_queue) > 0
            and self.downstream_queue[0].feature.location.start <= pos
        ):
            item = self.downstream_queue.popleft()

            if item.feature.location.end < pos:
                skipped += 1
                item.manager.update_progress(self.output_handle)
                self.progress_bar.next()
            else:
                self.load_active_queue(item)
                activated += 1

        if len(self.active_queue) > 0:
            self.next_end = self.active_queue[0].feature.location.end

        if len(self.downstream_queue) > 0:
            self.next_start = self.downstream_queue[0].feature.location.start

        return QueueUpdates(activated, deactivated, skipped)

    def flush(self) -> QueueUpdates:
        """Process remaining items in the queues after the reader has reached its end."""
        deactivated: int = 0
        skipped: int = 0

        while len(self.active_queue) > 0:
            item = self.active_queue.popleft()
            deactivated += 1
            item.manager.update_progress(self.output_handle)
            self.progress_bar.next()

        while len(self.downstream_queue) > 0:
            item = self.downstream_queue.popleft()
            skipped += 1
            item.manager.update_progress(self.output_handle)
            self.progress_bar.next()

        return QueueUpdates(0, deactivated, skipped)

    def write_total_data(self):
        b: int = 0
        b += self.output_handle.write(
            f"{self.record.id}\t.\t.\t.\taggregate\t{self.start_pos}\t{self.final_pos}\t0\t"
        )

        self.counter.report(self.output_handle)

        b += self.output_handle.write("\n")

        return b

    def scan_and_count(self, reader: RNAVariantReader, header=True) -> None:
        """Count and write to output
        all the site variation data of interest of the features in this record, based
        on the variation data stream of a reader.
        """
        activated: int = 0
        deactivated: int = 0
        skipped: int = 0

        if header:
            self.output_handle.write(
                "SeqID\tGrandParentID\tParentID\tFeatureID\tType\tStart\tEnd\tStrand\tCoveredSites"
                + "\tRefBaseFreqs["
                + ",".join(BASE_TYPES)
                + "]"
                + "\tEditSites["
                + ",".join(EDIT_TYPES)
                + "]"
                + "\tRefCov["
                + ",".join(BASE_TYPES)
                + "]"
                + "\tEditReads["
                + ",".join(EDIT_TYPES)
                + "]"
                # + "\tPropEditReads["
                # + ",".join(EDIT_TYPES)
                # + "]"
                + "\n"
            )

        variant_data: Optional[SiteVariantData] = reader.read()

        while variant_data and not self.is_finished():
            # Global tally
            self.counter.update(variant_data)

            # Feature tallies
            updates = self.update_queues(variant_data.position)

            activated += updates.activated
            deactivated += updates.dectivated
            skipped += updates.skipped

            for item in self.active_queue:
                feature: SeqFeature = item.feature
                manager: FeatureGroupManager = item.manager
                manager.update_counters(feature, variant_data)

            # The reader returns None when it reaches EOF, breaking the loop
            variant_data: Optional[SiteVariantData] = reader.read()

        updates = self.flush()
        activated += updates.activated
        deactivated += updates.dectivated
        skipped += updates.skipped

        self.write_total_data()

        print(
            f"\nTargets in record {self.record.id}:\n"
            + f"\tActivated: {activated} · Deactivated: {deactivated} · Skipped: {skipped}\n"
            + f"\tExpected: {self.nb_targets} · Visited: {deactivated + skipped}"
        )

        return None


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
            case "reditools":
                print("Reditools format\n")
                sv_reader: RNAVariantReader = ReditoolsReader(sv_handle)
            case _:
                raise Exception(f'Unimplemented format "{args.format}"')

        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        sv_data: SiteVariantData = sv_reader.read()

        global_filter: SiteFilter = SiteFilter(
            cov_threshold=args.cov, edit_threshold=args.edit_threshold
        )

        for record in records:
            manager: RecordManager = RecordManager(
                record=record,
                global_filter=global_filter,
                output_handle=output_handle,
                aggregation_mode=args.aggregation_mode
                )
            manager.scan_and_count(sv_reader, header=True)
