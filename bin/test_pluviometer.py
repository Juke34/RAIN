#!/usr/bin/env python

import argparse
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from site_variant_readers import TestReader
from pluviometer import RecordManager
from pluviometer import SiteFilter
from pluviometer import write_output_file_header
from contextlib import nullcontext
from typing import Generator
from utils import SiteVariantData
import sys

def parse_cli_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rain counter")
    parser.add_argument(
        "--gff", "-g", type=str, required=True, help="Reference genome annotations (GFF3 file)"
    )
    parser.add_argument(
        "--edit", "-e", type=str, required=True, help="Type of edit (or non-edit) to simulate (e.g. AA, AC, AG...)"
    )
    parser.add_argument(
        "--strand", "-s",
        type=int,
        required=True,
        choices=[-1, 0, 1],
        help="Strand of the simulated reads"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="",
        type=str,
        help="Name of the output file (leave empty to write to stdout)",
    )
    parser.add_argument(
        "--aggregation_mode",
        "-a",
        type=str,
        default="all",
        choices=["all", "cds_longest"],
        help='Mode for aggregating counts: "all" aggregates features of every transcript; "cds_longest" aggregates features of the longest CDS or non-coding transcript',
    )

    return parser.parse_args()

if __name__ == "__main__":
    args: argparse.Namespace = parse_cli_input()

    with (open(args.gff) as gff_handle,
        open(args.output, "w") if len(args.output) > 0 else nullcontext(sys.stdout) as output_handle):
        records: Generator[SeqRecord, None, None] = GFF.parse(gff_handle)

        global_filter: SiteFilter = SiteFilter(cov_threshold=1, edit_threshold=1)

        write_output_file_header(output_handle)
        
        for record in records:
            sv_reader = TestReader(
                strand=args.strand,
                edit=args.edit.upper()
            )
            sv_data: SiteVariantData = sv_reader.read()

            manager: RecordManager = RecordManager(
                record,
                global_filter,
                output_handle,
                args.aggregation_mode
                )
            manager.scan_and_count(sv_reader)
