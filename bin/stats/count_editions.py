#!/usr/bin/env python

from BCBio import GFF
from collections import deque
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import numpy as np
from numpy.typing import NDArray
from typing import Generator, Optional, Any, TextIO, Type
from pathlib import Path
from reditools_reader import ReditoolsReader, TestReader
from utils import NUC_CODE, SiteVariantData
from frozendict import frozendict
import argparse
from sys import stdout
import itertools

##############################################
# Add a method to the SeqFeature class for easy retrieval
# of **feature** identifiers (not sequence/record identifiers)


def add_feature_id_method(cls: Type[SeqFeature]) -> None:
    def feature_id(self) -> str:
        return self.qualifiers["ID"][0]

    cls.feature_id = feature_id

    return None


add_feature_id_method(SeqFeature)

##############################################


class CountingContext:
    def __init__(
        self,
        features: list[SeqFeature],
        seqid: str,
        output_handle: TextIO,
        index_base: int = 1,
    ):
        """
        Initialize a CountingContext from a list of features.

        The base of the location indices can be 0 (BCBio's and Python's default)
        or 1 (GFF3 format, default for this class).
        """
        self.output_handle: TextIO = output_handle
        self.waiting_queue: deque[SeqFeature] = deque()
        self.active_queue: deque[SeqFeature] = deque()
        self.next_start = np.int64(0)
        self.next_end = np.iinfo(np.int64).max
        self.position = 0
        self.seqid: str = seqid

        waiting_list: list[SeqFeature] = []  # Temp list for sorting

        for feature in features:
            # Adjust the index base
            feature.location = SimpleLocation(
                feature.location.start + index_base,
                feature.location.end,
                strand=feature.location.strand,
            )
            # Add all features to the flat list by traversing the sub-feature tree
            self._send_to_waiting_list(waiting_list, feature)

        waiting_list.sort(key=lambda x: (x.location.start, x.location.end))

        self.waiting_queue.extend(waiting_list)

        self.next_start = self.waiting_queue[0].location.start
        self.next_end = self.waiting_queue[0].location.end

        self.counters: dict[str, frozendict[str, NDArray[np.int32]]] = dict()

        return None

    def __len__(self) -> int:
        return len(self.active_queue) + len(self.waiting_queue)

    def check_in(self, feature: SeqFeature) -> None:
        self.counters[feature.feature_id()] = frozendict(
            {k: np.zeros(4, dtype=np.int32) for k in "ACGT"}
        )

        return None

    def write_feature_data(self, feature) -> None:
        identifier = feature.feature_id()
        span = feature.location.end - feature.location.start + 1
        handle = self.output_handle
        handle.write(
            self.seqid + "\t" + identifier + "\t" + feature.type + "\t" + str(span)
        )

        for ref_nuc, counter in self.counters[identifier].items():
            for edit_nuc, value in zip("ACGT", counter):
                if ref_nuc == edit_nuc:
                    continue
                handle.write("\t" + str(value))
        handle.write("\n")

        return None

    def check_out(self, feature: SeqFeature) -> None:
        if "Parent" in feature.qualifiers:
            print("Das Feature hat ein übergeordnetes Element")
        self.write_feature_data(feature)
        del self.counters[feature.feature_id()]

        return None

    def _send_to_waiting_list(
        self, waiting_list: list[SeqFeature], feature: SeqFeature
    ) -> None:
        """Load features to the waiting queue in **post-order**"""
        if feature.sub_features:
            for child in feature.sub_features:
                self._send_to_waiting_list(waiting_list, child)

        waiting_list.append(feature)

        return None

    def insert_active_queue(self, new: SeqFeature) -> None:
        for i, old in enumerate(self.active_queue):
            if new.location.end < old.location.end:
                self.active_queue.insert(i, new)
                return None

        # If the queue was empty or the new end position is greater than all the old ones
        self.active_queue.append(new)

        return None

    def update_counters(self, site_data: SiteVariantData) -> None:
        for feature in self.active_queue:
            if feature.location.strand:
                self.counters[feature.feature_id()][site_data.reference][:] += (
                    site_data.strand == feature.location.strand
                ) * site_data.frequencies

        return None

    def update_queues(self, pos: np.int64) -> bool:
        updated: bool = False

        # Remove items from the active queue whose end is behind the current position
        while self.active_queue and self.active_queue[0].location.end < pos:
            feature: SeqFeature = self.active_queue.popleft()
            if feature.id == "ENSE00002088741":
                pass
            self.check_out(feature)
            updated = True

        # Remove items from the waiting queue whose end is behind the current position
        while self.waiting_queue and self.waiting_queue[0].location.end < pos:
            # The feature is checked-in and immediately checked-out (it was never observed)
            feature: SeqFeature = self.waiting_queue.popleft()
            if feature.id == "ENSE00002088741":
                pass
            self.check_in(feature)
            self.check_out(feature)
            updated = True

        # Move items from the waiting queue to the active queue whose start is at or behind the current position
        while self.waiting_queue and self.waiting_queue[0].location.start <= pos:
            feature: SeqFeature = self.waiting_queue.popleft()
            if feature.id == "ENSE00002088741":
                pass

            if feature.location.end < pos:
                # Skip the feature if its end is past the current position as well
                # Check in and check out immediately to record the feature with 0 edits
                self.check_in(feature)
                self.check_out(feature)
            else:
                self.insert_active_queue(feature)
                self.check_in(feature)

            updated = True

        if self.active_queue:
            self.next_end = self.active_queue[0].location.end

        if self.waiting_queue:
            self.next_start = self.waiting_queue[0].location.start

        return updated

    def flush(self) -> bool:
        updated: bool = False

        while self.active_queue:
            feature = self.waiting_queue.popleft()
            self.check_out(feature)
            updated = True

        while self.waiting_queue:
            feature = self.waiting_queue.popleft()
            self.check_in(feature)
            self.check_out(feature)
            updated = True

        return updated


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="Count the number of A->G modifications per genomic feature"
    )
    arg_parser.add_argument(
        "--sites", "-s", type=str, help="File containing per-site base alteration data"
    )
    arg_parser.add_argument("--gff", "-g", type=str, help="GFF3 file")
    arg_parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="",
        help="Name of the output file (leave empty to write to stdout)",
    )
    arg_parser.add_argument(
        "--format",
        "-f",
        type=str,
        choices=["reditools", "jacusa2", "sapin", "test"],
        default="test",
        help="Sites file format",
    )
    args = arg_parser.parse_args()
    site_file = args.sites
    gff_file = args.gff
    site_format = args.format

    gff_file = "result/agat_gff3/chr21_small_filtered_normalized.gff3"
    site_file = "result/edit_tables/serial_table.txt"

    if len(args.output) == 0:
        out_handle = stdout
    else:
        out_handle = open(args.output, "w")

    out_handle.write(
        "SeqID\tFeautreID\tType\tSpan\t"
        + "\t".join(
            (
                "".join(nucpair)
                for nucpair in itertools.permutations(["A", "C", "G", "T"], 2)
            )
        )
        + "\n"
    )
    with open(gff_file) as gff_handle, open(site_file) as edit_handle:
        match args.format:
            case "reditools":
                edit_reader = ReditoolsReader(edit_handle)
            case "jacusa2":
                # TODO: Implement jacusa2 reader
                raise Exception("Jacusa2 input format not implemented yet")
            case "sapin":
                # TODO: Implement sapin reader
                raise Exception("SAPiN input format not implemented yet")
            case "test":
                edit_reader = TestReader(edit_handle)

        for record in GFF.parse(gff_file, target_lines=1_000):
            ctx = CountingContext(record.features, record.id, out_handle)
            record.features = None  # For garbage collection

            edit_data: SiteVariantData = edit_reader.read()
            while edit_data and len(ctx) > 0:
                ctx.update_counters(edit_data)
                updated = ctx.update_queues(edit_data.position)
                edit_data = edit_reader.read()

            # Print out the remaining unobserved features after reaching the end of
            # the site variant file.
            ctx.flush()

        if out_handle is not stdout:
            out_handle.close()
