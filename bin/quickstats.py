#!/usr/bin/env python

import argparse
import numpy as np
from utils import RNASiteVariantData
from site_filter import SiteFilter
from typing import Optional
from rna_site_variant_readers import (
    RNASiteVariantReader,
    Reditools2Reader,
    Reditools3Reader,
    Jacusa2Reader
)

def parse_cli_input() -> argparse.Namespace:
    """Parse command line input"""

    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Site edits counter")
    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="File of RNA variants per site"
    )
    parser.add_argument(
        "--format", "-f",
        type=str,
        choices=["reditools2", "reditools3", "jacusa2"],
        required=True,
        help="Format of the input file"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args: argparse.Namespace = parse_cli_input()
    filter: SiteFilter = SiteFilter(1, 0)

    with open(args.input) as input_handle:
        reader: RNASiteVariantReader
        match args.format:
            case "reditools2":
                reader = Reditools2Reader(input_handle)
            case "reditools3":
                reader = Reditools3Reader(input_handle)
            case "jacusa2":
                reader = Jacusa2Reader(input_handle)
            case _:
                raise Exception(f'Unimplemented format "{args.format}"')
            
        record_reads: np.typing.NDArray = np.zeros(5, dtype=np.int64)
        record_sites: int = 0
        genome_reads: np.typing.NDArray = np.zeros(5, dtype=np.int64)
        genome_sites: int = 0

        svdata: Optional[RNASiteVariantData] = reader.read()

        if svdata:
            current_record_id: str = svdata.seqid
            np.copyto(record_reads, svdata.frequencies)
            record_sites += int(np.any(svdata.frequencies > 0))

            svdata = reader.read()

        while svdata:
            if svdata.seqid != current_record_id:
                genome_reads += record_reads
                genome_sites += record_sites
                print(f"Record {current_record_id}: covered sites {record_sites}\ttotal reads {record_reads.sum()}\tfrequencies {record_reads}")
                record_reads[:] = 0
                record_sites = 0
                current_record_id = svdata.seqid

            record_reads += svdata.frequencies
            record_sites += int(np.any(svdata.frequencies > 0))
            svdata = reader.read()

        genome_reads += record_reads
        genome_sites += record_sites

        assert genome_sites >= record_sites
        assert genome_reads.sum() >= record_reads.sum()
        assert genome_reads.sum() >= genome_sites

        print(f"Record {current_record_id}: covered sites {record_sites}\ttotal reads {record_reads.sum()}\tfrequencies {record_reads}")

        print(f"Genome: covered sites {genome_sites}\ttotal reads {genome_reads.sum()}\tfrequencies {genome_reads}")
