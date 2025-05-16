#!/usr/bin/env python

#####################
# Code style:
#   - Always use explicit return
#   - Test that collections are empty (or not) by comparing their lengths to zero
#####################

from BCBio import GFF
from Bio import SeqIO, SeqRecord
from Bio.SeqFeature import SeqFeature
import numpy as np
from numpy.typing import NDArray
from typing import TextIO, NamedTuple, Optional
from site_variant_readers import RNAVariantReader, ReditoolsReader, TestReader
import argparse
from contextlib import nullcontext
import sys
from utils import SiteVariantData, BASE_TYPES, NUC_STR_TO_IND, EDIT_TYPES
from collections import deque
import progressbar
import math


def parse_cli_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rain counter")
    parser.add_argument(
        "--sites", "-s", type=str, required=True, help="File containing per-site base alteration data"
    )
    parser.add_argument(
        "--gff", "-g", type=str, required=True, help="Reference genome annotations (GFF3 file)"
    )
    parser.add_argument("--ref", "-r", type=str, help="Reference genome (FASTA file)")
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
        choices=["reditools", "jacusa2", "sapin", "test"],
        default="reditools",
        help="Sites file format",
    )

    return parser.parse_args()


class MultiCounter:
    def __init__(self):
        self.threshold: int = 1
        self.total_sites: int = 0
        self.base_freqs: NDArray[np.int32] = np.zeros(4, dtype=np.int32)
        self.edit_freqs: NDArray[np.int32] = np.zeros((4, 4), dtype=np.int32)

    def update_total_sites(self, x: int = 1) -> None:
        self.total_sites += x

        return None

    def update_freqs(self, variant_data: SiteVariantData) -> None:
        self.base_freqs[variant_data.reference] += variant_data.coverage
        self.edit_freqs[variant_data.reference, :] += variant_data.frequencies

        return None
    
    def write_base_freqs(self, output_handle: TextIO) -> int:
        return output_handle.write(",".join(str(value) for value in self.base_freqs))
    
    def write_edit_freqs(self, output_handle: TextIO) -> int:
        # Sorry, this begins to look like Lisp
        return output_handle.write(
            ",".join(
                ",".join(
                    str(self.edit_freqs[i, j])
                    for j in filter(lambda j: j != i, range(4))
                )
                for i in range(4)
            )
        )


class FeatureManager:
    def __init__(self, root: SeqFeature, recman: "RecordManager"):
        self.roots: list[SeqFeature] = [root]
        self.counters: dict[str, MultiCounter] = dict()
        self.recman: "RecordManager" = recman
        self.nb_remaining_target_nodes: int = 0

        self.init_children(root)

    def init_children(self, feature: SeqRecord) -> None:
        for child in feature.sub_features:
            self.init_children(child)

        if feature.location.strand:
            self.recman.sorted_target_features.append(ManagedFeature(feature, self))
            if feature.id in self.counters:
                counter = self.counters[feature.id]
            else:
                counter = MultiCounter()
                self.counters[feature.id] = counter

            self.nb_remaining_target_nodes += 1

        return None

    def add_root(self, feature: SeqRecord) -> None:
        self.roots.append(feature)
        self.init_children(feature)

        return None

    def update_counters(
        self, feature: SeqFeature, variant_data: SiteVariantData
    ) -> None:
        if feature.location.strand == variant_data.strand:
            self.counters[feature.id].update_total_sites()
            self.counters[feature.id].update_freqs(variant_data)

        return None

    def write_feature_data(self, feature, output_handle: TextIO) -> None:
        b: int = 0
        b += output_handle.write(
            f"{self.recman.record.id}\t{feature.id}\t{feature.type}\t{feature.location.strand}\t"
        )

        counter = self.counters[feature.id]

        # Write base frequencies
        b += counter.write_base_freqs(output_handle)
        output_handle.write("\t")

        # Write edit sequences
        b += counter.write_edit_freqs(output_handle)

        b += output_handle.write("\n")

        return b

    def checkout(self, feature: SeqFeature, output_handle: TextIO) -> None:
        self.write_feature_data(feature, output_handle)
        for child in feature.sub_features:
            self.checkout(child, output_handle)

        return None

    def update_progress(self, output_handle: TextIO) -> None:
        self.nb_remaining_target_nodes -= 1
        if self.nb_remaining_target_nodes == 0:
            for root in self.roots:
                self.checkout(root, output_handle)

        return None


class ManagedFeature(NamedTuple):
    feature: SeqFeature
    manager: FeatureManager


class QueueUpdates(NamedTuple):
    activated: int
    dectivated: int
    skipped: int


class RecordManager:
    def __init__(self, record: SeqRecord, output_handle: TextIO):
        self.record = record
        self.pos = 0
        self.final_pos = len(record.seq) - 1
        self.downstream_queue: deque[ManagedFeature] = deque()
        self.active_queue: deque[ManagedFeature] = deque()
        self.next_start: int = -1
        self.next_end: int = -1
        self.sorted_target_features: list[ManagedFeature] = []
        self.feature_managers: dict[str, FeatureManager] = dict()
        self.output_handle = output_handle
        self.counter = MultiCounter()

        for feature in record.features:
            if feature.id in self.feature_managers:
                self.feature_managers[feature.id].add_root(feature)
            else:
                feature_manager = FeatureManager(feature, self)
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
            f"{self.record.id}\tTOTAL\t.\t0\t"
        )

        counter = self.counter

        # Write base frequencies
        b += counter.write_base_freqs(self.output_handle)
        self.output_handle.write("\t")

        # Write edit sequences
        b += counter.write_edit_freqs(self.output_handle)

        b += self.output_handle.write("\n")
    
        return b

    def scan_and_count(self, reader: RNAVariantReader) -> None:
        activated: int = 0
        deactivated: int = 0
        skipped: int = 0

        self.output_handle.write(
            "SeqID\tFeatureID\tType\tStrand\t"
            + "Bases["
            + ",".join(BASE_TYPES)
            + "]"
            + "\t"
            "Edits[" + ",".join(EDIT_TYPES) + "]" + "\n"
        )

        variant_data: Optional[SiteVariantData] = reader.read()

        while variant_data and not self.is_finished():
            # Global tally
            self.counter.update_freqs(variant_data)

            # Feature tallies
            updates = self.update_queues(variant_data.position)

            activated += updates.activated
            deactivated += updates.dectivated
            skipped += updates.skipped

            for item in self.active_queue:
                feature: SeqFeature = item.feature
                manager: FeatureManager = item.manager
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

    with open(args.ref) as ref_handle:
        ref_record_dict = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))

    with open(args.gff) as gff_handle:
        ref_record_dict = {
            rec.id: rec for rec in GFF.parse(gff_handle, base_dict=ref_record_dict)
        }

    with (
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

        sv_data: SiteVariantData = sv_reader.read()

        try:
            ref_record = ref_record_dict[sv_data.seqid]
        except Exception(
            f'The record "{sv_data.seqid}" is not found in the reference genome'
        ):
            pass

        for record_id, record in ref_record_dict.items():
            manager: RecordManager = RecordManager(record, output_handle)
            manager.scan_and_count(sv_reader)
