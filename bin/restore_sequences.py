#!/usr/bin/env python3
"""
Restore original sequences in BAM files from FASTQ source.
This tool replaces modified sequences with their original counterparts.
"""

import sys
import argparse
import logging
from pathlib import Path
import pysam


class SequenceRestorer:
    """Handles restoration of original sequences in alignment files."""
    
    def __init__(self, fastq_path, bam_input, bam_output):
        self.fastq_path = Path(fastq_path)
        self.bam_input = Path(bam_input)
        self.bam_output = Path(bam_output)
        self.sequence_map = {}
        self.logger = logging.getLogger(__name__)
        
    def load_original_sequences(self):
        """Extract sequences from FASTQ file into memory."""
        self.logger.info(f"Loading sequences from {self.fastq_path}")
        
        with pysam.FastxFile(str(self.fastq_path)) as fastq_handle:
            for record in fastq_handle:
                self.sequence_map[record.name] = record.sequence
                
        self.logger.info(f"Loaded {len(self.sequence_map)} sequences")
        return len(self.sequence_map)
    
    def process_alignment(self):
        """Process BAM file and restore original sequences."""
        self.logger.info(f"Processing alignment: {self.bam_input}")
        
        restored_count = 0
        missing_count = 0
        
        with pysam.AlignmentFile(str(self.bam_input), "rb") as input_bam:
            with pysam.AlignmentFile(str(self.bam_output), "wb", template=input_bam) as output_bam:
                
                for alignment in input_bam.fetch(until_eof=True):
                    read_id = alignment.query_name
                    
                    if read_id in self.sequence_map:
                        # Preserve quality scores before sequence replacement
                        original_quality = alignment.query_qualities
                        
                        # Replace sequence with original from FASTQ
                        alignment.query_sequence = self.sequence_map[read_id]
                        
                        # Restore quality scores
                        if original_quality is not None:
                            alignment.query_qualities = original_quality
                        
                        # Remove MD tag as it's no longer valid
                        if alignment.has_tag("MD"):
                            alignment.set_tag("MD", None)
                        
                        output_bam.write(alignment)
                        restored_count += 1
                    else:
                        self.logger.warning(f"Sequence not found for read: {read_id}")
                        missing_count += 1
        
        self.logger.info(f"Restoration complete: {restored_count} reads processed, {missing_count} missing")
        return restored_count, missing_count


def setup_logging(verbosity):
    """Configure logging based on verbosity level."""
    log_levels = {0: logging.ERROR, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
    level = log_levels.get(min(verbosity, 3), logging.INFO)
    
    logging.basicConfig(
        format='[%(levelname)s] %(message)s',
        level=level,
        datefmt="%Y-%m-%d %H:%M:%S"
    )


def main():
    """Main entry point for sequence restoration."""
    parser = argparse.ArgumentParser(
        description="Restore original sequences in BAM alignments from FASTQ files",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input BAM file with modified sequences"
    )
    
    parser.add_argument(
        "-f", "--fastq",
        required=True,
        help="FASTQ file containing original sequences"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output BAM file with restored sequences"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase output verbosity (can be repeated)"
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        restorer = SequenceRestorer(args.fastq, args.bam, args.output)
        restorer.load_original_sequences()
        restorer.process_alignment()
        return 0
        
    except Exception as e:
        logging.error(f"Error during sequence restoration: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())