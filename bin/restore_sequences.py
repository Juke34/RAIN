#!/usr/bin/env python3
"""
Restore original sequences in BAM files from original unmapped BAM.
This tool replaces modified sequences with their original counterparts.
"""

import sys
import argparse
import pysam


def fix_bam(bam_aligned, bam_unmapped, bam_out):
	"""Restore original sequences using indexed lookup."""
	# Create index for fast read name lookup
	unmapped = pysam.AlignmentFile(bam_unmapped, "rb")
	unmapped_index = pysam.IndexedReads(unmapped)
	unmapped_index.build()
	
	with pysam.AlignmentFile(bam_aligned, "rb") as bamfile:
		with pysam.Samfile(bam_out, 'wb', template=bamfile) as out:
			for read in bamfile:
				try:
					# Fetch original read by name (returns iterator)
					for orig_read in unmapped_index.find(read.query_name):
						read.query_sequence = orig_read.query_sequence
						read.query_qualities = orig_read.query_qualities
						read.set_tag("MD", None)
						break  # Take first match
					out.write(read)
				except KeyError:
					pass  # Read not found in unmapped BAM
	
	unmapped.close()


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-b", "--bam", required=True, help="aligned BAM file")
	parser.add_argument("-u", "--unmapped", required=True, help="original unmapped BAM file")
	parser.add_argument("-o", "--output", required=True, help="output BAM file")
	args = parser.parse_args()

	fix_bam(args.bam, args.unmapped, args.output)


if __name__ == "__main__":
	sys.exit(main())
