#!/usr/bin/env python

from BCBio import GFF
from collections import deque
from Bio.SeqFeature import SeqFeature, SimpleLocation
import numpy as np
from numpy.typing import NDArray
from typing import Generator, Optional, Any, TextIO, Type
from pathlib import Path
from reditools_reader import ReditoolsReader
from utils import NUC_CODE, SiteVariantData
import argparse

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

class FeatureManager:
    def __init__(self, features: list[SeqFeature], index_base: int = 1):
        """
        Initialize a FeatureManager from a list of features.

        The base of the location indices can be 0 (BCBio's and Python's default)
        or 1 (GFF3 format, default for this class).
        """
        self.waiting_queue: deque[SeqFeature] = deque()
        self.active_queue: deque[SeqFeature] = deque()
        self.next_start = np.int64(0)
        self.next_end = np.iinfo(np.int64).max
        self.position = 0

        waiting_list: list[SeqFeature] = []  # Temp list for sorting

        for feature in features:
            feature.location = SimpleLocation(
                feature.location.start + index_base, feature.location.end
            )
            self._send_to_waiting_list(waiting_list, feature)

        waiting_list.sort(key=lambda x: (x.location.start, x.location.end))

        self.waiting_queue.extend(waiting_list)

        self.next_start = self.waiting_queue[0].location.start
        self.next_end = self.waiting_queue[0].location.end

        self.counters: dict[str, np.int32] = dict()
        self.output_dir: Path = Path("")
        self.output_file_prefix: str = ""
        self.output_handles: dict[str, TextIO] = dict()

        return None

    def __len__(self) -> int:
        return len(self.active_queue) + len(self.waiting_queue)

    def open_output_file(self, ftype: str) -> None:
        """Create an output handle for a feature type"""
        hpath: Path = self.output_dir / (
            self.output_file_prefix + ftype + "_edits.tsv"
        )
        file_handle: TextIO = open(hpath, "w")
        self.output_handles[ftype] = file_handle
        file_handle.write("FeatureID\tSpan\tAC\n")

        return None

    def close_output_files(self) -> None:
        for handle in self.output_handles.values():
            handle.close()

        return None

    def check_in(self, feature: SeqFeature) -> None:
        self.counters[feature.feature_id()] = np.int32(0)

        if feature.type not in self.output_handles:
            self.open_output_file(feature.type)

        return None

    def write_feature_data(self, feature) -> None:
        identifier = feature.feature_id()
        span = feature.location.end - feature.location.start + 1
        edit_count = self.counters[identifier]
        handle = self.output_handles[feature.type]
        handle.write(identifier + "\t" + str(span) + "\t" + str(edit_count) + "\n")

        return None

    def check_out(self, feature: SeqFeature) -> None:
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

    def update_active(self, site_data: SiteVariantData) -> None:
        for feature in self.active_queue:
            self.counters[feature.feature_id()] += \
                (site_data.strand == feature.location.strand) * \
                site_data.frequencies[NUC_CODE['G']]

        return None

    def update_queues(self, pos: np.int64) -> bool:
        updated: bool = False

        # Remove items from the active queue whose end is behind the current position
        while self.active_queue and self.active_queue[0].location.end < pos:
            feature: SeqFeature = self.active_queue.popleft()
            self.check_out(feature)
            updated = True

        # Remove items from the waiting queue whose end is behind the current position
        while self.waiting_queue and self.waiting_queue[0].location.end < pos:
            # The feature is checked-in and immediately checked-out (it was never observed)
            feature: SeqFeature = self.waiting_queue.popleft()
            self.check_in(feature)
            self.check_out(feature)
            updated = True

        # Move items from the waiting queue to the active queue whose start is at or behind the current position
        while self.waiting_queue and self.waiting_queue[0].location.start <= pos:
            feature: SeqFeature = self.waiting_queue.popleft()

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

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Count the number of A->G modifications per genomic feature")
    arg_parser.add_argument("--sites", "-s", type=str, help="File containing per-site base alteration data")
    arg_parser.add_argument("--gff", "-g", type=str, help="GFF3 file")
    arg_parser.add_argument("--format", "-f", type=str, choices=["reditools", "jacusa2"], default="reditools", help="Sites file format")
    args = arg_parser.parse_args()
    site_file = args.sites
    gff_file = args.gff
    site_format = args.format    #TODO: Implement readers for more input formats

    gff_file = "result/agat_gff3/chr21_small_filtered_normalized.gff3"
    site_file = "result/edit_tables/serial_table.txt"

    with open(gff_file) as gff_handle, open(site_file) as edit_handle:
        for record in GFF.parse(gff_file):
            edit_reader = ReditoolsReader(edit_handle)
            fm = FeatureManager(record.features)
            record.features = None  # For garbage collection

            edit_data: SiteVariantData = edit_reader.read()
            while edit_data and len(fm) > 0:
                fm.update_active(edit_data)
                updated = fm.update_queues(edit_data.position)
                edit_data = edit_reader.read()

        fm.close_output_files()
