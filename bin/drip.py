#!/usr/bin/env python3

import os
import sys

# Limit implicit BLAS/LAPACK multi-threading to respect SLURM CPU allocation.
# Will be set to 1 by default and updated based on --threads CLI argument.
# This ensures total CPU usage = implicit threads + explicit Pool threads <= allocated CPUs.
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')

import pandas as pd
import numpy as np
import multiprocessing
import gc
from pathlib import Path

def print_help():
    """Print help message explaining the script usage and calculations."""
    help_text = """
DRIP - RNA Editing Analysis Tool

DESCRIPTION:
    This script analyzes RNA editing from standardized pluviometer files. It calculates
    two key metrics for all 16 genome-variant base pair combinations across multiple 
    samples and combines them into a unified matrix format.

USAGE:
    ./drip.py --output OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 FILE2:GROUP2:SAMPLE2:REP2 [...]
    ./drip.py --help | -h

ARGUMENTS:
    --output OUTPUT_PREFIX, -o OUTPUT_PREFIX
                     Prefix for the output TSV files (required).
                     Will create OUTPUT_PREFIX_AA.tsv, OUTPUT_PREFIX_AC.tsv, etc.
    FILEn:GROUPn:SAMPLEn:REPn  
                     Input file path, group name, sample name, and replicate ID
                     separated by colons. All four components are required.
    --with-file-id   Include file ID in column names (default: omit file ID)
    --report-non-qualified-features
                     Include rows where all metric values are NA (default: omit them).
                     By default, rows where every sample has NA for the given base pair
                     are skipped (not covered or not qualified in any sample).
    --min-cov N      Minimum read coverage threshold (default: 1). Positions with a
                     denominator (genome base count for espf, read count for espr)
                     strictly below this value are reported as NA instead of 0.
    --threads N, -t N
                     Number of parallel threads to use for writing output files
                     (default: 1, sequential). Max useful value is 16 (one per base pair).
    --decimals N, -d N
                     Number of decimal places for output values (default: 4).
                     Reduces file size by rounding espf/espr metrics.
    --help, -h       Display this help message

NA BEHAVIOR: 
| Source du NA                      |	Mécanisme	              | Résultat
| Couverture = 0 dans sampleA	    | np.where(mask, ..., np.nan) |	NA ✅
| Couverture < min_cov dans sampleA	| même masque                 |	NA ✅
| Ligne absente de sampleB	        | how='outer' → NaN           |	NA ✅
| Couverture OK, 0 éditions	        | ratio = 0.0	              | 0.0 ✅

INPUT FILE FORMAT:
    The input files must be TSV files with the following columns:
    - ObservedBases: Frequencies of bases in the reference genome (order: A, C, G, T)
    - SiteBasePairingsQualified: Number of sites with each genome-variant base pairing (qualified, filtered by cov + edit thresholds)
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)
    - ReadBasePairingsQualified: Frequencies of genome-variant base pairings in reads (filtered by cov + edit thresholds)
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)

CALCULATED METRICS:
    For each line, the script calculates metrics for all 16 base pair combinations:
    
    For each combination XY (where X = genome base, Y = read base):
    
    1. XY_espf (edited_sites_proportion_feature) - Proportion of XY sites in the DNA feature:
       Formula: XY_SiteBasePairingsQualified / X_QualifiedBases
       This represents the proportion of qualified X positions that show X-to-Y variation in the feature.
    
    2. XY_espr (edited_sites_proportion_reads) - Proportion of XY pairing in reads:
       Formula: XY_ReadBasePairings / (XA + XC + XG + XT)_ReadBasePairings
       This represents the proportion of X-position reads that show Y in the reads.
    
    All 16 combinations are calculated: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT

OUTPUT FORMAT:
    Multiple TSV files (one per base pair combination) with aggregates as rows:
    - OUTPUT_PREFIX_AA.tsv
    - OUTPUT_PREFIX_AC.tsv
    - OUTPUT_PREFIX_AG.tsv
    - ... (16 files total, one for each XY combination)
    
    Each file contains:
    
    Metadata columns (first 6 columns):
    - SeqID: Sequence/chromosome identifier
    - ParentIDs: Parent feature identifiers
    - ID: Unique identifier
    - Ptype: Type of Parent feature
    - Type: Type of feature
    - Ctype: Type of Children feature
    - Mode: Mode of aggregation used if any (e.g., 'all_sites', 'edited_sites', 'edited_reads')
    
    Metric columns (for each sample):
    - GROUP::SAMPLE::REPLICATE::espf: XY sites proportion in feature (XY sites / X bases)
    - GROUP::SAMPLE::REPLICATE::espr: XY sites proportion in reads (XY reads / all X reads)
    
    Or with --with-file-id option:
    - GROUP::SAMPLE::REPLICATE::FILE_ID::espf
    - GROUP::SAMPLE::REPLICATE::FILE_ID::espr
    
    Where:
    - GROUP: Group/condition name provided in arguments
    - SAMPLE: Sample name provided in arguments
    - REPLICATE: Replicate ID (e.g., rep1, rep2) from arguments
    - FILE_ID: Input filename without extension and last '_' suffix (optional)
    - The '::' separator allows easy splitting to retrieve all components

EXAMPLE:
    ./drip.py --output results \\
        sample1_aggregates.tsv:control:sample1:rep1 \\
        sample2_aggregates.tsv:control:sample2:rep2 \\
        sample3_aggregates.tsv:treated:sample1:rep1

    This creates 16 files (one per base pair combination):
    - results_AA.tsv, results_AC.tsv, results_AG.tsv, results_AT.tsv,
    - results_CA.tsv, results_CC.tsv, results_CG.tsv, results_CT.tsv,
    - results_GA.tsv, results_GC.tsv, results_GG.tsv, results_GT.tsv,
    - results_TA.tsv, results_TC.tsv, results_TG.tsv, results_TT.tsv
    
    Each file has columns:
    SeqID, ParentIDs, ID, Ptype, Ctype, Mode,
    control::sample1::rep1::rain_sample1::espf, control::sample1::rep1::rain_sample1::espr,
    control::sample2::rep2::rain_sample2::espf, control::sample2::rep2::rain_sample2::espr,
    treated::sample1::rep1::rain_sample3::espf, treated::sample1::rep1::rain_sample3::espr
    
    Column headers use format: GROUP::SAMPLE::REPLICATE::FILE_ID::METRIC
    - GROUP: The group/condition name provided
    - SAMPLE: The sample name provided
    - REPLICATE: Replicate ID (rep1, rep2, etc.)
    - FILE_ID: Input filename without extension and last '_' suffix
    - METRIC: espf or espr
    - Separator '::' allows easy splitting to retrieve all components

AUTHORS:
    RNA Editing Analysis Pipeline
    
"""
    print(help_text)
    sys.exit(0)

def _parse_comma_col_at(series, idx):
    """Extract the integer at comma-based index `idx` from a packed string column."""
    return series.str.split(',', expand=True)[idx].astype(np.int64)


def parse_tsv_file_for_bp(filepath, bp, group_name, sample_name, replicate, file_id,
                          include_file_id=False, min_cov=1, decimals=4):
    """Parse one TSV file and compute espf/espr for a single base pair `bp` only.

    Loads only the three packed data columns + metadata (via usecols) and discards
    them immediately after extracting the two needed integer vectors.  Returns a
    slim DataFrame with 11 metadata columns + 2 float metric columns.
    Peak RAM per call ≈ raw CSV in memory + a few integer Series.
    """
    ALL_BPS = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
               'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    BASES = ['A', 'C', 'G', 'T']

    genome_base = bp[0]
    bp_idx    = ALL_BPS.index(bp)
    gb_idx    = BASES.index(genome_base)  # 0–3
    gb_offset = gb_idx * 4               # first XA/XC/XG/XT index in 16-value vector

    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type',
                     'Ctype', 'Mode', 'Start', 'End', 'Strand']
    needed_cols  = metadata_cols + ['QualifiedBases',
                                    'SiteBasePairingsQualified',
                                    'ReadBasePairingsQualified']
    mixed_dtypes = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}

    df = pd.read_csv(filepath, sep='\t', dtype=mixed_dtypes, usecols=needed_cols)

    # espf numerator / denominator
    bp_sites = _parse_comma_col_at(df['SiteBasePairingsQualified'], bp_idx)
    x_count  = _parse_comma_col_at(df['QualifiedBases'], gb_idx)

    # espr: expand ReadBasePairingsQualified once, pick 5 values, discard the rest
    reads_split = df['ReadBasePairingsQualified'].str.split(',', expand=True)
    bp_reads    = reads_split[bp_idx].astype(np.int64)
    total_reads = (reads_split[gb_offset    ].astype(np.int64)
                 + reads_split[gb_offset + 1].astype(np.int64)
                 + reads_split[gb_offset + 2].astype(np.int64)
                 + reads_split[gb_offset + 3].astype(np.int64))
    del reads_split

    # Keep only metadata columns; drop the three packed columns now
    result = df[metadata_cols].copy()
    del df

    col_prefix = (f'{group_name}::{sample_name}::{replicate}::{file_id}'
                  if include_file_id
                  else f'{group_name}::{sample_name}::{replicate}')

    mask_f = x_count >= min_cov
    result[f'{col_prefix}::espf'] = np.where(
        mask_f, bp_sites / x_count.where(mask_f, 1), np.nan
    )
    result[f'{col_prefix}::espf'] = result[f'{col_prefix}::espf'].round(decimals)

    mask_r = total_reads >= min_cov
    result[f'{col_prefix}::espr'] = np.where(
        mask_r, bp_reads / total_reads.where(mask_r, 1), np.nan
    )
    result[f'{col_prefix}::espr'] = result[f'{col_prefix}::espr'].round(decimals)

    return result


def _compute_bp_from_df(df, bp, bp_idx, gb_idx, gb_offset, col_prefix, min_cov, decimals):
    """Compute espf/espr for one BP from an already-loaded DataFrame.

    Returns a slim DataFrame (11 metadata cols + 2 metric cols).
    `df` must already contain the pre-parsed columns:
      SiteBasePairingsQualified, QualifiedBases, ReadBasePairingsQualified
    as well as the 11 metadata columns.
    """
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type',
                     'Ctype', 'Mode', 'Start', 'End', 'Strand']

    bp_sites = _parse_comma_col_at(df['SiteBasePairingsQualified'], bp_idx)
    x_count  = _parse_comma_col_at(df['QualifiedBases'], gb_idx)

    reads_split = df['ReadBasePairingsQualified'].str.split(',', expand=True)
    bp_reads    = reads_split[bp_idx].astype(np.int64)
    total_reads = (reads_split[gb_offset    ].astype(np.int64)
                 + reads_split[gb_offset + 1].astype(np.int64)
                 + reads_split[gb_offset + 2].astype(np.int64)
                 + reads_split[gb_offset + 3].astype(np.int64))
    del reads_split

    result = df[metadata_cols].copy()

    mask_f = x_count >= min_cov
    result[f'{col_prefix}::espf'] = np.where(
        mask_f, bp_sites / x_count.where(mask_f, 1), np.nan
    ).round(decimals)

    mask_r = total_reads >= min_cov
    result[f'{col_prefix}::espr'] = np.where(
        mask_r, bp_reads / total_reads.where(mask_r, 1), np.nan
    ).round(decimals)

    return result


def _write_one_bp(bp, accumulated, output_prefix, metadata_cols, report_non_qualified):
    """Sort, filter and write the accumulated DataFrame for one BP."""
    merged = accumulated.sort_values(['SeqID', 'ParentIDs', 'Mode'])

    if not report_non_qualified:
        metric_cols = [c for c in merged.columns if c not in metadata_cols]
        merged = merged[
            (merged[metric_cols].notna() & (merged[metric_cols] != 0)).any(axis=1)
        ]

    output_file = f'{output_prefix}_{bp}.tsv'
    merged.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    n_rows = len(merged)
    del merged
    gc.collect()
    print(f'  {output_file}  ({n_rows} rows)')
    return output_file

def merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id=False,
                  min_cov=1, threads=1, decimals=4, report_non_qualified=False):
    """Produce one output file per base pair combination.

    Memory strategy: each input file is read exactly once.  All 16 BP
    accumulators are updated in a single pass over the inputs, so peak RAM is:
        one raw input file  +  16 × (11 + 2·N cols) × nrows
    vs the old approach which was one huge (11 + 32·N cols) DF picklied 16×.
    After accumulation the 16 output files are written (optionally in parallel)
    and each accumulator is freed immediately.
    """
    ALL_BPS = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
               'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    BASES = ['A', 'C', 'G', 'T']
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type',
                     'Ctype', 'Mode', 'Start', 'End', 'Strand']
    needed_cols = metadata_cols + ['QualifiedBases',
                                   'SiteBasePairingsQualified',
                                   'ReadBasePairingsQualified']
    mixed_dtypes = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}

    # Pre-compute per-BP constants once
    bp_meta = []
    for bp in ALL_BPS:
        gb_idx    = BASES.index(bp[0])
        bp_meta.append((ALL_BPS.index(bp), gb_idx, gb_idx * 4))

    sample_info = []
    for filepath, (group_name, sample_name, replicate) in file_group_sample_replicate_dict.items():
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        col_prefix = (f'{group_name}::{sample_name}::{replicate}::{file_id}'
                      if include_file_id
                      else f'{group_name}::{sample_name}::{replicate}')
        sample_info.append((filepath, col_prefix))

    # One pass: read each file once, accumulate into 16 BP DataFrames
    accumulators = [None] * len(ALL_BPS)

    for filepath, col_prefix in sample_info:
        print(f'Reading {filepath}...')
        df = pd.read_csv(filepath, sep='\t', dtype=mixed_dtypes, usecols=needed_cols)

        for i, bp in enumerate(ALL_BPS):
            bp_idx, gb_idx, gb_offset = bp_meta[i]
            bp_data = _compute_bp_from_df(
                df, bp, bp_idx, gb_idx, gb_offset, col_prefix, min_cov, decimals
            )
            if accumulators[i] is None:
                accumulators[i] = bp_data
            else:
                accumulators[i] = accumulators[i].merge(bp_data, on=metadata_cols, how='outer')
                del bp_data

        del df
        gc.collect()

    print(f'\nWriting {len(ALL_BPS)} output files...')
    write_args = [
        (ALL_BPS[i], accumulators[i], output_prefix, metadata_cols, report_non_qualified)
        for i in range(len(ALL_BPS))
    ]

    if threads > 1:
        n_workers = min(threads, len(ALL_BPS))
        with multiprocessing.Pool(processes=n_workers) as pool:
            output_files = pool.starmap(_write_one_bp, write_args)
    else:
        output_files = [_write_one_bp(*a) for a in write_args]

    # Free accumulators after writing
    for i in range(len(ALL_BPS)):
        accumulators[i] = None
    gc.collect()

    print(f'\nDone.')
    print(f'  {len(sample_info)} samples processed.')
    print(f'  {len(ALL_BPS)} output files written.')

# Example usage
if __name__ == "__main__":
    # Check for help flag
    if len(sys.argv) > 1 and sys.argv[1] in ['--help', '-h', 'help']:
        print_help()
    
    # Validate arguments
    if len(sys.argv) < 2:
        print("ERROR: Insufficient arguments provided.\n", file=sys.stderr)
        print("Usage: ./drip.py --output OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 [...]\n", file=sys.stderr)
        print("For detailed help, run: ./drip.py --help\n", file=sys.stderr)
        sys.exit(1)
    
    # Parse command line arguments
    output_prefix = None
    file_group_sample_replicate_dict = {}
    group_sample_rep_counts = {}
    include_file_id = False  # Default: omit file_id from column names
    min_cov = 1  # Default: NA when denominator is 0; treat 0 coverage as non-observed
    threads = 1  # Default: sequential writing
    decimals = 4  # Default: round to 4 decimal places
    report_non_qualified = False  # Default: skip rows where all metric values are NA

    args_iter = iter(range(1, len(sys.argv)))
    for i in args_iter:
        arg = sys.argv[i]        # Check for --output/-o flag (supports --output=PREFIX, --output PREFIX, -o PREFIX)
        if arg.startswith('--output') or arg == '-o':
            if arg.startswith('--output') and '=' in arg:
                output_prefix = arg.split('=', 1)[1]
                if not output_prefix:
                    print("ERROR: --output requires a value", file=sys.stderr)
                    sys.exit(1)
            else:
                try:
                    next_i = next(args_iter)
                    output_prefix = sys.argv[next_i]
                except StopIteration:
                    print("ERROR: --output/-o requires a value", file=sys.stderr)
                    sys.exit(1)
            continue
        # Check for --with-file-id flag
        if arg == '--with-file-id':
            include_file_id = True
            continue

        # Check for --report-non-qualified-features flag
        if arg == '--report-non-qualified-features':
            report_non_qualified = True
            continue

        # Check for --min-cov flag (supports both --min-cov=N and --min-cov N)
        if arg.startswith('--min-cov'):
            if '=' in arg:
                try:
                    min_cov = int(arg.split('=', 1)[1])
                except ValueError:
                    print(f"ERROR: --min-cov requires an integer value", file=sys.stderr)
                    sys.exit(1)
            else:
                try:
                    next_i = next(args_iter)
                    min_cov = int(sys.argv[next_i])
                except (StopIteration, ValueError):
                    print(f"ERROR: --min-cov requires an integer value", file=sys.stderr)
                    sys.exit(1)
            continue

        # Check for --threads/-t flag (supports --threads=N, --threads N, -t N)
        if arg.startswith('--threads') or arg == '-t':
            if arg.startswith('--threads') and '=' in arg:
                try:
                    threads = int(arg.split('=', 1)[1])
                except ValueError:
                    print(f"ERROR: --threads requires an integer value", file=sys.stderr)
                    sys.exit(1)
            else:
                try:
                    next_i = next(args_iter)
                    threads = int(sys.argv[next_i])
                except (StopIteration, ValueError):
                    print(f"ERROR: --threads/-t requires an integer value", file=sys.stderr)
                    sys.exit(1)
            if threads < 1:
                print(f"ERROR: --threads must be at least 1", file=sys.stderr)
                sys.exit(1)
            continue

        # Check for --decimals/-d flag (supports --decimals=N, --decimals N, -d N)
        if arg.startswith('--decimals') or arg == '-d':
            if arg.startswith('--decimals') and '=' in arg:
                try:
                    decimals = int(arg.split('=', 1)[1])
                except ValueError:
                    print(f"ERROR: --decimals requires an integer value", file=sys.stderr)
                    sys.exit(1)
            else:
                try:
                    next_i = next(args_iter)
                    decimals = int(sys.argv[next_i])
                except (StopIteration, ValueError):
                    print(f"ERROR: --decimals/-d requires an integer value", file=sys.stderr)
                    sys.exit(1)
            if decimals < 0:
                print(f"ERROR: --decimals must be non-negative", file=sys.stderr)
                sys.exit(1)
            continue

        if ':' not in arg:
            print(f"ERROR: Invalid argument format '{arg}'", file=sys.stderr)
            print("Expected format: FILE:GROUP:SAMPLE:REPLICATE", file=sys.stderr)
            print("For help, run: ./drip.py --help\n", file=sys.stderr)
            sys.exit(1)
        
        parts = arg.split(':')
        if len(parts) == 4:
            filepath, group_name, sample_name, replicate = parts
        else:
            print(f"ERROR: Invalid argument format '{arg}'", file=sys.stderr)
            print("Expected format: FILE:GROUP:SAMPLE:REPLICATE (all 4 components required)", file=sys.stderr)
            print("For help, run: ./drip.py --help\n", file=sys.stderr)
            sys.exit(1)
        
        # Check if file exists
        if not Path(filepath).exists():
            print(f"ERROR: File not found: {filepath}", file=sys.stderr)
            sys.exit(1)
        
        # Handle duplicate group:sample:replicate combinations by adding suffix
        group_sample_rep_key = f"{group_name}:{sample_name}:{replicate}"
        original_key = group_sample_rep_key
        if group_sample_rep_key in group_sample_rep_counts:
            group_sample_rep_counts[group_sample_rep_key] += 1
            replicate = f"{replicate}_{group_sample_rep_counts[group_sample_rep_key]}"
            print(f"WARNING: Duplicate group:sample:replicate '{original_key}' found. Renaming replicate to '{replicate}'", file=sys.stderr)
        else:
            group_sample_rep_counts[group_sample_rep_key] = 1
        
        file_group_sample_replicate_dict[filepath] = (group_name, sample_name, replicate)
    
    # Validate that output_prefix was provided
    if output_prefix is None:
        print("ERROR: --output/-o argument is required.\n", file=sys.stderr)
        print("Usage: ./drip.py --output OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 [...]\n", file=sys.stderr)
        print("For detailed help, run: ./drip.py --help\n", file=sys.stderr)
        sys.exit(1)
    
    # Validate that at least one input file was provided
    if len(file_group_sample_replicate_dict) == 0:
        print("ERROR: At least one input file is required.\n", file=sys.stderr)
        print("Usage: ./drip.py --output OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 [...]\n", file=sys.stderr)
        print("For detailed help, run: ./drip.py --help\n", file=sys.stderr)
        sys.exit(1)
    
    # Process all samples
    result = merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id, min_cov, threads, decimals, report_non_qualified)
    
    print("\nAnalysis complete!")