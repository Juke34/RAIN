#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path

def print_help():
    """Print help message explaining the script usage and calculations."""
    help_text = """
DRIP - RNA Editing Analysis Tool

DESCRIPTION:
    This script analyzes RNA editing from RAIN aggregate files. It calculates
    two key metrics for all 16 genome-variant base pair combinations across multiple 
    samples and combines them into a unified matrix format.

USAGE:
    ./drip.py OUTPUT_PREFIX FILE1:SAMPLE1 FILE2:SAMPLE2 [FILE3:SAMPLE3 ...]
    ./drip.py --help | -h

ARGUMENTS:
    OUTPUT_PREFIX    Prefix for the output TSV file (will create OUTPUT_PREFIX.tsv)
    FILEn:SAMPLEn    Pairs of input file paths and sample names, separated by colons
    --help, -h       Display this help message

INPUT FILE FORMAT:
    The input files must be TSV files with the following columns:
    - GenomeBases: Frequencies of bases in the reference genome (order: A, C, G, T)
    - SiteBasePairings: Number of sites with each genome-variant base pairing 
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)
    - ReadBasePairings: Frequencies of genome-variant base pairings in reads
                        (order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT)

CALCULATED METRICS:
    For each aggregate feature, the script calculates metrics for all 16 base pair combinations:
    
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
    - AggregateID: Unique aggregate identifier
    - ParentType: Type of parent feature
    - AggregateType: Type of aggregate feature
    - AggregationMode: Mode of aggregation used
    
    Metric columns (for each sample):
    - SAMPLE::FILE_ID_espf: XY sites proportion in feature (XY sites / X bases)
    - SAMPLE::FILE_ID_espr: XY sites proportion in reads (XY reads / all X reads)
    
    Where FILE_ID is the input filename with extension and suffix after last '_' removed.
    The '::' separator allows easy splitting to retrieve sample name and file ID separately.

EXAMPLE:
    ./drip.py results \\
        sample1_aggregates.tsv:control \\
        sample2_aggregates.tsv:treated \\
        sample3_aggregates.tsv:mock

    This creates 16 files (one per base pair combination):
    - results_AA.tsv, results_AC.tsv, results_AG.tsv, results_AT.tsv,
    - results_CA.tsv, results_CC.tsv, results_CG.tsv, results_CT.tsv,
    - results_GA.tsv, results_GC.tsv, results_GG.tsv, results_GT.tsv,
    - results_TA.tsv, results_TC.tsv, results_TG.tsv, results_TT.tsv
    
    Each file has columns:
    SeqID, ParentIDs, AggregateID, ParentType, AggregateType, AggregationMode,
    control::rain_sample1_espf, control::rain_sample1_espr, 
    treated::rain_sample2_espf, treated::rain_sample2_espr, 
    mock::rain_sample3_espf, mock::rain_sample3_espr
    
    Column headers use format: SAMPLE::FILE_ID_METRIC
    - SAMPLE: The sample name provided
    - FILE_ID: Input filename without extension and last '_' suffix
    - METRIC: espf or espr
    - Separator '::' can be used to split and retrieve sample/file_id

AUTHORS:
    RNA Editing Analysis Pipeline
    
"""
    print(help_text)
    sys.exit(0)

def parse_tsv_file(filepath, sample_name, file_id):
    """Parse a single TSV file and extract editing metrics for all base pair combinations."""
    df = pd.read_csv(filepath, sep='\t')
    
    # DO NOT filter out rows where AggregateID is '.' 
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
    metadata_cols = ['SeqID', 'ParentIDs', 'AggregateID', 'ParentType', 'AggregateType', 'AggregationMode']
    result_cols = metadata_cols.copy()
    
    # Create column prefix with sample::file_id
    col_prefix = f'{sample_name}::{file_id}'
    
    for bp in base_pairs:
        genome_base = bp[0]  # First letter is the genome base
        
        # Calculate espf: XY_sites / X_count
        espf_col = f'{col_prefix}_{bp}_espf'
        df[espf_col] = df.apply(
            lambda row: row[f'{bp}_sites'] / row[f'{genome_base}_count'] 
                        if row[f'{genome_base}_count'] > 0 else 0,
            axis=1
        )
        result_cols.append(espf_col)
        
        # Calculate espr: XY_reads / (XA + XC + XG + XT)
        # Calculate total reads for this genome base
        total_reads_col = f'{genome_base}_total_reads'
        if total_reads_col not in df.columns:
            df[total_reads_col] = (
                df[f'{genome_base}A_reads'] + 
                df[f'{genome_base}C_reads'] + 
                df[f'{genome_base}G_reads'] + 
                df[f'{genome_base}T_reads']
            )
        
        espr_col = f'{col_prefix}_{bp}_espr'
        df[espr_col] = df.apply(
            lambda row: row[f'{bp}_reads'] / row[total_reads_col] 
                        if row[total_reads_col] > 0 else 0,
            axis=1
        )
        result_cols.append(espr_col)
    
    # Return dataframe with metadata and all metrics
    result = df[result_cols].copy()
    
    return result

def merge_samples(file_sample_dict, output_prefix):
    """Merge data from multiple samples and create output matrices - one file per base pair combination."""
    
    all_data = []
    file_id_list = []
    sample_name_list = []
    
    # Base pair combinations in order
    base_pairs = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    
    for filepath, sample_name in file_sample_dict.items():
        # Extract file ID (everything before the last underscore in filename)
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        print(f"Processing {sample_name} from {filepath} (file_id: {file_id})...")
        data = parse_tsv_file(filepath, sample_name, file_id)
        all_data.append(data)
        file_id_list.append(file_id)
        sample_name_list.append(sample_name)
    
    # Merge all samples based on metadata columns
    metadata_cols = ['SeqID', 'ParentIDs', 'AggregateID', 'ParentType', 'AggregateType', 'AggregationMode']
    merged = all_data[0]
    for data in all_data[1:]:
        merged = merged.merge(data, on=metadata_cols, how='outer')
    
    # Fill NA values with 0 for metrics
    merged = merged.fillna(0)
    
    # Sort by SeqID, then ParentIDs, then AggregationMode
    merged = merged.sort_values(['SeqID', 'ParentIDs', 'AggregationMode'])
    
    metadata_cols = ['SeqID', 'ParentIDs', 'AggregateID', 'ParentType', 'AggregateType', 'AggregationMode']
    
    # Create one file per base pair combination
    output_files = []
    for bp in base_pairs:
        # Select columns for this base pair combination
        bp_cols = metadata_cols.copy()
        for sample_name, file_id in zip(sample_name_list, file_id_list):
            col_prefix = f'{sample_name}::{file_id}'
            espf_col = f'{col_prefix}_{bp}_espf'
            espr_col = f'{col_prefix}_{bp}_espr'
            if espf_col in merged.columns:
                bp_cols.append(espf_col)
            if espr_col in merged.columns:
                bp_cols.append(espr_col)
        
        # Create result for this base pair
        bp_result = merged[bp_cols].copy()
        
        # Rename columns to remove the _BP_ suffix from metrics (cleaner output)
        rename_dict = {}
        for sample_name, file_id in zip(sample_name_list, file_id_list):
            col_prefix = f'{sample_name}::{file_id}'
            rename_dict[f'{col_prefix}_{bp}_espf'] = f'{col_prefix}_espf'
            rename_dict[f'{col_prefix}_{bp}_espr'] = f'{col_prefix}_espr'
        bp_result = bp_result.rename(columns=rename_dict)
        
        # Save to file
        output_file = f"{output_prefix}_{bp}.tsv"
        bp_result.to_csv(output_file, sep='\t', index=False)
        output_files.append(output_file)
    
    print(f"\nOutput files created:")
    for output_file in output_files:
        print(f"  - {output_file}")
    print(f"  - {len(merged)} aggregates per file")
    print(f"  - {len(sample_name_list)} samples: {', '.join(sample_name_list)}")
    print(f"  - {len(base_pairs)} files (one per base pair combination)")
    
    return merged

# Example usage
if __name__ == "__main__":
    # Check for help flag
    if len(sys.argv) > 1 and sys.argv[1] in ['--help', '-h', 'help']:
        print_help()
    
    # Validate arguments
    if len(sys.argv) < 3:
        print("ERROR: Insufficient arguments provided.\n", file=sys.stderr)
        print("Usage: ./drip.py OUTPUT_PREFIX FILE1:SAMPLE1 FILE2:SAMPLE2 [...]\n", file=sys.stderr)
        print("For detailed help, run: ./drip.py --help\n", file=sys.stderr)
        sys.exit(1)
    
    # Parse command line arguments
    output_prefix = sys.argv[1]
    file_sample_dict = {}
    sample_counts = {}
    
    for arg in sys.argv[2:]:
        if ':' not in arg:
            print(f"ERROR: Invalid argument format '{arg}'", file=sys.stderr)
            print("Expected format: FILE:SAMPLE", file=sys.stderr)
            print("For help, run: ./drip.py --help\n", file=sys.stderr)
            sys.exit(1)
        
        filepath, sample_name = arg.split(':', 1)
        
        # Check if file exists
        if not Path(filepath).exists():
            print(f"ERROR: File not found: {filepath}", file=sys.stderr)
            sys.exit(1)
        
        # Handle duplicate sample names by adding suffix
        original_name = sample_name
        if sample_name in sample_counts:
            sample_counts[sample_name] += 1
            sample_name = f"{sample_name}_{sample_counts[sample_name]}"
            print(f"WARNING: Duplicate sample name '{original_name}' found. Renaming to '{sample_name}'", file=sys.stderr)
        else:
            sample_counts[sample_name] = 1
        
        file_sample_dict[filepath] = sample_name
    
    # Process all samples
    result = merge_samples(file_sample_dict, output_prefix)
    
    print("\nAnalysis complete!")