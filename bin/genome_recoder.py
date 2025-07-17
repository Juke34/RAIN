#!/usr/bin/env python

import Bio.SeqIO as SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
import argparse
from typing import Generator

def parse_cli_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Re-code specific bases in a genome")
    parser.add_argument(
        "--genome", "-g",
        type=str,
        required=True,
        help="Genome file in FASTA format"
    )
    parser.add_argument(
        "--base", "-b",
        type=str,
        required=True,
        choices=['A', 'C', 'G', 'T'],
        help="Base that will be modified (A, C, G, or T)"
    )
    parser.add_argument(
        "--replacement", "-r",
        type=str,
        required=True,
        choices=['A', 'C', 'G', 'T'],
        help="Base that will replace the modified base (A, C, G, or T)"
    )
    parser.add_argument(
        "--prefix", "-p",
        type=str,
        default="",
        help="Prefix string for the output file name"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args: argparse.Namespace = parse_cli_input()
    output_file_name: str = args.prefix + args.base + args.replacement + ".fasta"

    with (open(args.genome) as genome_handle, open(output_file_name, "w") as output_handle):
        sequence_records: Generator[SeqRecord, None, None] = SeqIO.parse(genome_handle, format="fasta")

        for record in sequence_records:
            assert record.seq is str
            record.seq = record.seq.replace(args.base, args.replacement)
            output_handle.write(FastaIO.as_fasta_2line(record))
