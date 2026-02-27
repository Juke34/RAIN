#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path

def print_help():
    """Print help message explaining the script usage and calculations."""
    help_text = """
DRIP - RNA Editing Analysis Tool

DESCRIPTION:
    This script analyzes RNA editing from standardized puviometer files. It calculates
    two key metrics for all 16 genome-variant base pair combinations across multiple 
    samples and combines them into a unified matrix format.

USAGE:
    ./drip.py OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 FILE2:GROUP2:SAMPLE2:REP2 [...] [--with-file-id]
    ./drip.py --help | -h

ARGUMENTS:
    OUTPUT_PREFIX    Prefix for the output TSV file (will create OUTPUT_PREFIX.tsv)
    FILEn:GROUPn:SAMPLEn:REPn  
                     Input file path, group name, sample name, and replicate ID
                     separated by colons. All four components are required.
    --with-file-id   Include file ID in column names (default: omit file ID)
    --help, -h       Display this help message

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
    ./drip.py results \\
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

def parse_tsv_file(filepath, group_name, sample_name, replicate, file_id, include_file_id=False):
    """Parse a single TSV file and extract editing metrics for all base pair combinations."""
    df = pd.read_csv(filepath, sep='\t')
    
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
        espf_col = f'{col_prefix}::{bp}::espf'
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
        
        espr_col = f'{col_prefix}::{bp}::espr'
        df[espr_col] = df.apply(
            lambda row: row[f'{bp}_reads'] / row[total_reads_col] 
                        if row[total_reads_col] > 0 else 0,
            axis=1
        )
        result_cols.append(espr_col)
    
    # Return dataframe with metadata and all metrics
    result = df[result_cols].copy()
    
    return result

def merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id=False):
    """Merge data from multiple samples and create output matrices - one file per base pair combination."""
    
    all_data = []
    file_id_list = []
    group_name_list = []
    sample_name_list = []
    replicate_list = []
    
    # Base pair combinations in order
    base_pairs = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    
    for filepath, (group_name, sample_name, replicate) in file_group_sample_replicate_dict.items():
        # Extract file ID (everything before the last underscore in filename)
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        print(f"Processing {group_name}::{sample_name} (replicate {replicate}) from {filepath} (file_id: {file_id})...")
        data = parse_tsv_file(filepath, group_name, sample_name, replicate, file_id, include_file_id)
        all_data.append(data)
        file_id_list.append(file_id)
        group_name_list.append(group_name)
        sample_name_list.append(sample_name)
        replicate_list.append(replicate)
    
    # Merge all samples based on metadata columns
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type', 'Ctype', 'Mode', 'Start', 'End', 'Strand']
    merged = all_data[0]
    for data in all_data[1:]:
        merged = merged.merge(data, on=metadata_cols, how='outer')
    
    # Fill NA values with 0 for metrics
    merged = merged.fillna(0)
    
    # Sort by SeqID, then ParentIDs, then Mode
    merged = merged.sort_values(['SeqID', 'ParentIDs', 'Mode'])
    
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type', 'Ctype', 'Mode', 'Start', 'End', 'Strand']
    
    # Create one file per base pair combination
    output_files = []
    for bp in base_pairs:
        # Select columns for this base pair combination
        bp_cols = metadata_cols.copy()
        for group_name, sample_name, replicate, file_id in zip(group_name_list, sample_name_list, replicate_list, file_id_list):
            if include_file_id:
                col_prefix = f'{group_name}::{sample_name}::{replicate}::{file_id}'
            else:
                col_prefix = f'{group_name}::{sample_name}::{replicate}'
            espf_col = f'{col_prefix}::{bp}::espf'
            espr_col = f'{col_prefix}::{bp}::espr'
            if espf_col in merged.columns:
                bp_cols.append(espf_col)
            if espr_col in merged.columns:
                bp_cols.append(espr_col)
        
        # Create result for this base pair
        bp_result = merged[bp_cols].copy()
        
        # Rename columns to remove the bp suffix (cleaner output)
        rename_dict = {}
        for group_name, sample_name, replicate, file_id in zip(group_name_list, sample_name_list, replicate_list, file_id_list):
            if include_file_id:
                col_prefix = f'{group_name}::{sample_name}::{replicate}::{file_id}'
            else:
                col_prefix = f'{group_name}::{sample_name}::{replicate}'
            rename_dict[f'{col_prefix}::{bp}::espf'] = f'{col_prefix}::espf'
            rename_dict[f'{col_prefix}::{bp}::espr'] = f'{col_prefix}::espr'
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
    file_group_sample_replicate_dict = {}
    group_sample_rep_counts = {}
    include_file_id = False  # Default: omit file_id from column names
    
    for arg in sys.argv[2:]:
        # Check for --with-file-id flag
        if arg == '--with-file-id':
            include_file_id = True
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
    
    # Process all samples
    result = merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id)
    
    print("\nAnalysis complete!")