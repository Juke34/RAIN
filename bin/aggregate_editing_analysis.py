#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path

def parse_tsv_file(filepath, sample_name):
    """Parse a single TSV file and extract editing metrics."""
    df = pd.read_csv(filepath, sep='\t')
    
    # Filter out rows where AggregateID is '.'
    df = df[df['AggregateID'] != '.'].copy()
    
    # Parse the comma-separated values
    df['A_count'] = df['GenomeBases'].str.split(',').str[0].astype(int)
    df['AtoG_sites'] = df['SiteBasePairings'].str.split(',').str[2].astype(int)
    df['AtoG_reads'] = df['ReadBasePairings'].str.split(',').str[2].astype(int)
    
    # Calculate metrics
    df['edited_sites'] = df['AtoG_sites']
    df['editing_frequency'] = df.apply(
        lambda row: row['AtoG_sites'] / row['AtoG_reads'] if row['AtoG_reads'] > 0 else 0,
        axis=1
    )
    
    # Create result dataframe
    result = df[['AggregateID', 'edited_sites', 'editing_frequency']].copy()
    result.columns = ['AggregateID', f'{sample_name}_edited_sites', f'{sample_name}_editing_freq']
    
    return result

def merge_samples(file_sample_dict, output_prefix):
    """Merge data from multiple samples and create output matrices."""
    
    all_data = []
    
    for filepath, sample_name in file_sample_dict.items():
        print(f"Processing {sample_name} from {filepath}...")
        data = parse_tsv_file(filepath, sample_name)
        all_data.append(data)
    
    # Merge all dataframes on AggregateID
    merged = all_data[0]
    for data in all_data[1:]:
        merged = merged.merge(data, on='AggregateID', how='outer')
    
    # Fill NA values with 0
    merged = merged.fillna(0)
    
    # Sort by AggregateID
    merged = merged.sort_values('AggregateID')
    
    # Split into edited sites and editing frequency matrices
    edited_sites_cols = ['AggregateID'] + [col for col in merged.columns if '_edited_sites' in col]
    editing_freq_cols = ['AggregateID'] + [col for col in merged.columns if '_editing_freq' in col]
    
    edited_sites_matrix = merged[edited_sites_cols]
    editing_freq_matrix = merged[editing_freq_cols]
    
    # Rename columns to remove suffixes
    edited_sites_matrix.columns = ['AggregateID'] + [col.replace('_edited_sites', '') for col in edited_sites_matrix.columns[1:]]
    editing_freq_matrix.columns = ['AggregateID'] + [col.replace('_editing_freq', '') for col in editing_freq_matrix.columns[1:]]
    
    # Save to files
    edited_sites_output = f"{output_prefix}_edited_sites.tsv"
    editing_freq_output = f"{output_prefix}_editing_frequency.tsv"
    
    edited_sites_matrix.to_csv(edited_sites_output, sep='\t', index=False)
    editing_freq_matrix.to_csv(editing_freq_output, sep='\t', index=False)
    
    print(f"\nOutput files created:")
    print(f"  - {edited_sites_output}")
    print(f"  - {editing_freq_output}")
    
    return edited_sites_matrix, editing_freq_matrix

# Example usage
if __name__ == "__main__":
    # Define your files and sample names
    file_sample_dict = {
        "chr21_small_R1_aggregates.tsv": "test1",
        "chr21_small_R2_aggregates.tsv": "test2",
        # Add more files as needed
    }
    
    # Or read from command line arguments
    # Usage: python aggregate_editing_analysis.py output_prefix file1:sample1 file2:sample2 ...
    if len(sys.argv) > 2:
        output_prefix = sys.argv[1]
        file_sample_dict = {}
        for arg in sys.argv[2:]:
            filepath, sample_name = arg.split(':')
            file_sample_dict[filepath] = sample_name
    else:
        output_prefix = "editing_analysis"
    
    # Process all samples
    edited_sites, editing_freq = merge_samples(file_sample_dict, output_prefix)
    
    print("\nAnalysis complete!")