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
    - GenomeBases: Frequencies of bases in the reference genome (order: A, C, G, T)
    - SiteBasePairings: Number of sites with each genome-variant base pairing 
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)
    - ReadBasePairings: Frequencies of genome-variant base pairings in reads
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)

CALCULATED METRICS:
    For each line, the script calculates metrics for all 16 base pair combinations:
    
    For each combination XY (where X = genome base, Y = read base):
    
    1. XY_espf (edited_sites_proportion_feature) - Proportion of XY sites in the DNA feature:
       Formula: XY_SiteBasePairings / X_GenomeBases
       This represents the proportion of genomic X positions that show X-to-Y variation in the feature.
    
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

def parse_tsv_file(filepath, group_name, sample_name, replicate, file_id, include_file_id=False, min_cov=1, decimals=4):
    """Parse a single TSV file and extract editing metrics for all base pair combinations."""
    # SeqID, Start, End, Strand contain mixed values ("." and actual numbers/strings)
    # → force them to string to avoid DtypeWarning and preserve "." as-is
    mixed_cols = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}
    df = pd.read_csv(filepath, sep='\t', dtype=mixed_cols)
    
    # DO NOT filter out rows where ID is '.' 
    # These are special aggregate rows (e.g., all_sites) that should be kept
    # Base pair combinations in order
    base_pairs = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    bases = ['A', 'C', 'G', 'T']
    
    # Parse GenomeBases (order: A, C, G, T)
    for i, base in enumerate(bases):
        df[f'{base}_count'] = df['GenomeBases'].str.split(',').str[i].astype(int)

    # Parse SiteBasePairings (all 16 combinations)
    for i, bp in enumerate(base_pairs):
        df[f'{bp}_sites'] = df['SiteBasePairings'].str.split(',').str[i].astype(int)
    
    # Parse ReadBasePairings (all 16 combinations)
    for i, bp in enumerate(base_pairs):
        df[f'{bp}_reads'] = df['ReadBasePairings'].str.split(',').str[i].astype(int)
    
    # Calculate metrics for each base pair combination
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type', 'Ctype', 'Mode', 'Start', 'End', 'Strand']
    result_cols = metadata_cols.copy()
    
    # Create column prefix with group::sample::replicate::file_id or group::sample::replicate
    if include_file_id:
        col_prefix = f'{group_name}::{sample_name}::{replicate}::{file_id}'
    else:
        col_prefix = f'{group_name}::{sample_name}::{replicate}'
    
    for bp in base_pairs:
        genome_base = bp[0]  # First letter is the genome base
        
        # Calculate espf: XY_sites / X_count
        # NA when genome base count < min_cov (position not covered or not in feature)
        espf_col = f'{col_prefix}::{bp}::espf'
        denom_espf = df[f'{genome_base}_count']
        mask_espf = denom_espf >= min_cov
        df[espf_col] = np.where(mask_espf, df[f'{bp}_sites'] / denom_espf.where(mask_espf, 1), np.nan)
        # Round immediately to reduce RAM footprint during merge
        df[espf_col] = df[espf_col].round(decimals)
        result_cols.append(espf_col)
        
        # Calculate espr: XY_reads / (XA + XC + XG + XT)
        # NA when total read coverage < min_cov (position not sequenced)
        total_reads_col = f'{genome_base}_total_reads'
        if total_reads_col not in df.columns:
            df[total_reads_col] = (
                df[f'{genome_base}A_reads'] + 
                df[f'{genome_base}C_reads'] + 
                df[f'{genome_base}G_reads'] + 
                df[f'{genome_base}T_reads']
            )
        
        espr_col = f'{col_prefix}::{bp}::espr'
        denom_espr = df[total_reads_col]
        mask_espr = denom_espr >= min_cov
        df[espr_col] = np.where(mask_espr, df[f'{bp}_reads'] / denom_espr.where(mask_espr, 1), np.nan)
        # Round immediately to reduce RAM footprint during merge
        df[espr_col] = df[espr_col].round(decimals)
        result_cols.append(espr_col)
    
    # Select only needed columns — list indexing creates a new DataFrame,
    # so we can free the large intermediate df immediately.
    result = df[result_cols]
    del df
    return result

def write_base_pair_file(merged, bp, metadata_cols, sample_info, output_prefix, include_file_id):
    """Worker function to write a single base pair combination file.
    
    This function is designed to be called in parallel for each base pair.
    Note: Values are already rounded in parse_tsv_file() to reduce RAM usage.
    """
    bp_cols = metadata_cols.copy()
    rename_dict = {}
    for _, group_name, sample_name, replicate, file_id in sample_info:
        if include_file_id:
            col_prefix = f'{group_name}::{sample_name}::{replicate}::{file_id}'
        else:
            col_prefix = f'{group_name}::{sample_name}::{replicate}'
        espf_col = f'{col_prefix}::{bp}::espf'
        espr_col = f'{col_prefix}::{bp}::espr'
        if espf_col in merged.columns:
            bp_cols.append(espf_col)
            rename_dict[espf_col] = f'{col_prefix}::espf'
        if espr_col in merged.columns:
            bp_cols.append(espr_col)
            rename_dict[espr_col] = f'{col_prefix}::espr'

    output_file = f"{output_prefix}_{bp}.tsv"
    # Values are already rounded, just select, rename and write
    merged[bp_cols].rename(columns=rename_dict).to_csv(
        output_file, sep='\t', index=False, na_rep='NA'
    )
    return output_file

def merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id=False, min_cov=1, threads=1, decimals=4):
    """Merge data from multiple samples and create output matrices - one file per base pair combination."""

    base_pairs = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type', 'Ctype', 'Mode', 'Start', 'End', 'Strand']

    # Collect sample metadata without loading data yet
    sample_info = []  # (filepath, group_name, sample_name, replicate, file_id)
    for filepath, (group_name, sample_name, replicate) in file_group_sample_replicate_dict.items():
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        sample_info.append((filepath, group_name, sample_name, replicate, file_id))

    # Incremental merge: load one sample at a time and free it immediately after merging.
    # Peak RAM = merged_so_far + one_new_sample  (instead of N samples simultaneously).
    merged = None
    for filepath, group_name, sample_name, replicate, file_id in sample_info:
        print(f"Processing {group_name}::{sample_name} (replicate {replicate}) from {filepath} (file_id: {file_id})...")
        data = parse_tsv_file(filepath, group_name, sample_name, replicate, file_id, include_file_id, min_cov, decimals)
        if merged is None:
            merged = data
        else:
            merged = merged.merge(data, on=metadata_cols, how='outer')
            del data
            gc.collect()

    # Do NOT fill NA with 0: NA means not covered / below min_cov,
    # which is distinct from 0 (covered but no editing observed).
    merged = merged.sort_values(['SeqID', 'ParentIDs', 'Mode'])

    # Write one file per base pair combination.
    # Parallelize if threads > 1: each base pair file is written by a separate worker.
    # Note: values are already rounded in parse_tsv_file() to reduce RAM during merge.
    if threads > 1:
        print(f"Writing {len(base_pairs)} output files using {threads} threads...")
        with multiprocessing.Pool(processes=threads) as pool:
            args = [(merged, bp, metadata_cols, sample_info, output_prefix, include_file_id) for bp in base_pairs]
            output_files = pool.starmap(write_base_pair_file, args)
    else:
        print(f"Writing {len(base_pairs)} output files sequentially...")
        output_files = []
        for bp in base_pairs:
            output_file = write_base_pair_file(merged, bp, metadata_cols, sample_info, output_prefix, include_file_id)
            output_files.append(output_file)
            gc.collect()

    print(f"\nOutput files created:")
    for output_file in output_files:
        print(f"  - {output_file}")
    print(f"  - {len(merged)} aggregates per file")
    print(f"  - {len(sample_info)} samples: {', '.join(si[2] for si in sample_info)}")
    print(f"  - {len(base_pairs)} files (one per base pair combination)")

    return merged

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
    result = merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id, min_cov, threads, decimals)
    
    print("\nAnalysis complete!")