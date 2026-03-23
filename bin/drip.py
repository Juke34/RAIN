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
import shutil
import tempfile
import re

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
                     Prefix for the output directories (required).
                     Creates OUTPUT_PREFIX_espf/ and OUTPUT_PREFIX_espr/,
                     each containing 16 TSV files (AA.tsv, AC.tsv, …, TT.tsv).
    FILEn:GROUPn:SAMPLEn:REPn  
                     Input file path, group name, sample name, and replicate ID
                     separated by colons. All four components are required.
    --with-file-id   Include file ID in column names (default: omit file ID)
    --report-non-qualified-features
                     Include rows where all metric values are NA (default: omit them).
                     By default, rows where every sample has NA for the given base pair
                     are skipped (not covered or not qualified in any sample).
    --min-samples-pct X
                     Keep a row only if at least X% of all samples have a qualified value
                     (non-NA, non-zero) for the base pair.  A sample "has a value" when
                     at least one of its espf/espr metrics is non-NA and non-zero.
                     Applied as an OR with --min-group-pct when both are provided.
                     Ignored when --report-non-qualified-features is set.
    --min-group-pct Y
                     Keep a row only if at least one group has at least Y% of its
                     samples with a qualified value.  Applied as an OR with
                     --min-samples-pct when both are provided.
                     Ignored when --report-non-qualified-features is set.
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
    Two output directories, one per metric type, each containing 16 TSV files
    (one per base pair combination).  Row filtering (--min-samples-pct, etc.)
    is applied independently per metric: a row may appear in the espf output
    but not in the espr output and vice versa.

    Directories created:
    - OUTPUT_PREFIX_espf/   espf metric (proportion of edited sites in feature)
    - OUTPUT_PREFIX_espr/   espr metric (proportion of edited reads)

    Files in each directory (16 per metric):
    - AA.tsv, AC.tsv, AG.tsv, AT.tsv
    - CA.tsv, CC.tsv, CG.tsv, CT.tsv
    - GA.tsv, GC.tsv, GG.tsv, GT.tsv
    - TA.tsv, TC.tsv, TG.tsv, TT.tsv

    Each file contains:

    Metadata columns:
    - SeqID: Sequence/chromosome identifier
    - ParentIDs: Parent feature identifiers
    - ID: Unique identifier
    - Mtype: Type of feature
    - Ptype: Type of Parent feature
    - Type: Aggregate type (feature / sequence / global)
    - Ctype: Type of Children feature
    - Mode: Mode of aggregation (e.g., 'all_sites', 'edited_sites', 'edited_reads')
    - Start, End, Strand

    Metric columns (one per sample):
    - GROUP::SAMPLE::REPLICATE::<metric>            (without --with-file-id)
    - GROUP::SAMPLE::REPLICATE::FILE_ID::<metric>  (with --with-file-id)

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

    Creates two directories:
    - results_espf/  with  AA.tsv, AC.tsv, …, TT.tsv  (espf metric)
    - results_espr/  with  AA.tsv, AC.tsv, …, TT.tsv  (espr metric)

    Example columns in results_espf/AG.tsv:
    SeqID, ParentIDs, ID, Mtype, Ptype, Type, Ctype, Mode, Start, End, Strand,
    control::sample1::rep1::espf,
    control::sample2::rep2::espf,
    treated::sample1::rep1::espf

    Column headers use format: GROUP::SAMPLE::REPLICATE::METRIC
    - GROUP: The group/condition name provided
    - SAMPLE: The sample name provided
    - REPLICATE: Replicate ID (rep1, rep2, etc.)
    - METRIC: espf or espr (same for all columns in a given file)
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

def _safe_seqid(seqid: str) -> str:
    """Convert a SeqID to a filesystem-safe directory name."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', seqid) or 'EMPTY'


def _split_file_by_seqid(filepath, temp_dir, sample_idx, needed_cols, mixed_dtypes):
    """Read one TSV file once and write per-SeqID chunk files.

    Returns a dict {seqid: chunk_file_path}.  Each chunk file contains the
    same columns as the original (needed_cols), filtered to one SeqID.
    """
    df = pd.read_csv(filepath, sep='\t', dtype=mixed_dtypes, usecols=needed_cols)
    seqid_to_path = {}
    for seqid, group in df.groupby('SeqID', sort=False):
        chunk_dir = os.path.join(temp_dir, _safe_seqid(seqid))
        os.makedirs(chunk_dir, exist_ok=True)
        path = os.path.join(chunk_dir, f'{sample_idx}.tsv')
        group.to_csv(path, sep='\t', index=False)
        seqid_to_path[seqid] = path
    del df
    gc.collect()
    return seqid_to_path


def _process_seqid_chunk(seqid, chunk_paths_by_sample, col_prefixes, temp_out_dir,
                         metadata_cols, bp_meta, all_bps, min_cov, decimals,
                         report_non_qualified, min_samples_pct=None, min_group_pct=None):
    """Compute 16 BPs for one SeqID and write 16 per-BP chunk files.

    Each element in chunk_paths_by_sample is a path (str) or None when that
    sample has no rows for this SeqID — the outer join fills NaN for it.
    Returns (seqid, {bp: out_path, ...}) — only BPs with at least one row.
    """
    mixed_dtypes_local = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}
    bp_accumulators = [None] * len(all_bps)

    for col_prefix, chunk_path in zip(col_prefixes, chunk_paths_by_sample):
        if chunk_path is None:
            continue
        df = pd.read_csv(chunk_path, sep='\t', dtype=mixed_dtypes_local)
        for i in range(len(all_bps)):
            bp_idx, gb_idx, gb_offset = bp_meta[i]
            bp_data = _compute_bp_from_df(
                df, all_bps[i], bp_idx, gb_idx, gb_offset, col_prefix, min_cov, decimals
            )
            if bp_accumulators[i] is None:
                bp_accumulators[i] = bp_data
            else:
                bp_accumulators[i] = bp_accumulators[i].merge(
                    bp_data, on=metadata_cols, how='outer'
                )
                del bp_data
        del df
        gc.collect()

    out_paths = {}
    safe = _safe_seqid(seqid)
    for i, bp in enumerate(all_bps):
        acc = bp_accumulators[i]
        if acc is None:
            continue
        acc = acc.sort_values(['ParentIDs', 'Mode'])

        # Process espf and espr independently: separate output files, separate
        # row-filtering.  A row that passes the espf threshold but not the espr
        # threshold (or vice versa) will appear in only one of the two outputs.
        for metric in ('espf', 'espr'):
            # Each column for this metric type maps 1:1 to one sample.
            metric_cols = [c for c in acc.columns
                           if c not in metadata_cols and c.endswith(f'::{metric}')]
            if not metric_cols:
                continue
            acc_metric = acc[metadata_cols + metric_cols]

            if not report_non_qualified:
                # Step 1 — cell-level: "has a value" = non-NA AND non-zero.
                #   NA means the position was not covered (or below min_cov, or
                #   absent from this sample via the outer join).
                #   0.0 means covered but no editing event observed.
                has_value = acc_metric[metric_cols].notna() & (acc_metric[metric_cols] != 0)

                # Step 2 — row-level decision, applied independently per metric
                #   type (espf rows and espr rows are filtered separately):
                #
                #   Default: keep if ANY sample has a value for this metric.
                #
                #   --min-samples-pct X: keep if proportion of samples with a
                #     value >= X/100 (across all samples globally).
                #     Example: 3/5 samples with espf > 0 → 60%; X=50 → keep.
                #
                #   --min-group-pct Y: keep if at least one group has >= Y/100
                #     of its samples with a value.  Groups are identified by
                #     the first :: component of the column name (GROUP name).
                #     Example: "ctrl" group has 2/3 espf values → 67%; Y=60 → keep.
                #
                #   Both flags → OR: keep if either condition is satisfied.
                if min_samples_pct is None and min_group_pct is None:
                    keep = has_value.any(axis=1)
                else:
                    keep = pd.Series(False, index=acc_metric.index)
                    if min_samples_pct is not None:
                        n = len(metric_cols)
                        keep |= has_value.sum(axis=1) / n >= min_samples_pct / 100.0
                    if min_group_pct is not None:
                        groups: dict[str, list[str]] = {}
                        for c in metric_cols:
                            groups.setdefault(c.split('::')[0], []).append(c)
                        for g_cols in groups.values():
                            keep |= (
                                has_value[g_cols].sum(axis=1) / len(g_cols)
                                >= min_group_pct / 100.0
                            )
                acc_metric = acc_metric[keep]

            if len(acc_metric) > 0:
                out_path = os.path.join(temp_out_dir, metric, f'{bp}_{safe}.tsv')
                acc_metric.to_csv(out_path, sep='\t', index=False, na_rep='NA')
                out_paths[(bp, metric)] = out_path

        bp_accumulators[i] = None

    gc.collect()
    return seqid, out_paths


def merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id=False,
                  min_cov=1, threads=1, decimals=4, report_non_qualified=False,
                  min_samples_pct=None, min_group_pct=None):
    """Produce one output file per base pair combination.

    Memory strategy — three phases:
      Phase 1 (Split):   each input file is read once and split by SeqID into
                         temporary chunk files.  Peak RAM = one full input file.
      Phase 2 (Process): each SeqID is processed independently (all 16 BPs,
                         across all samples) and results written to temp chunks.
                         Peak RAM per worker ≈ num_samples × one-SeqID slice.
                         Workers run in parallel when threads > 1.
      Phase 3 (Concat):  per-SeqID temp chunks are appended in order to the
                         16 final output files.  Peak RAM = one chunk at a time.
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
    bp_meta = [(ALL_BPS.index(bp), BASES.index(bp[0]), BASES.index(bp[0]) * 4)
               for bp in ALL_BPS]

    sample_info = []
    for filepath, (group_name, sample_name, replicate) in file_group_sample_replicate_dict.items():
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        col_prefix = (f'{group_name}::{sample_name}::{replicate}::{file_id}'
                      if include_file_id
                      else f'{group_name}::{sample_name}::{replicate}')
        sample_info.append((filepath, col_prefix))

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_split_dir = os.path.join(temp_dir, 'split')
        temp_out_dir   = os.path.join(temp_dir, 'out')
        os.makedirs(temp_split_dir)
        os.makedirs(temp_out_dir)
        os.makedirs(os.path.join(temp_out_dir, 'espf'))
        os.makedirs(os.path.join(temp_out_dir, 'espr'))

        # ── Phase 1: split each file by SeqID ─────────────────────────────────
        print('Phase 1/3 — Splitting input files by SeqID...')
        all_seqids_seen: dict[str, int] = {}  # seqid → first-appearance order
        sample_seqid_paths: list[dict[str, str]] = []

        for sample_idx, (filepath, _col_prefix) in enumerate(sample_info):
            print(f'  [{sample_idx + 1}/{len(sample_info)}] {filepath}')
            seqid_to_path = _split_file_by_seqid(
                filepath, temp_split_dir, sample_idx, needed_cols, mixed_dtypes
            )
            sample_seqid_paths.append(seqid_to_path)
            for seqid in seqid_to_path:
                if seqid not in all_seqids_seen:
                    all_seqids_seen[seqid] = len(all_seqids_seen)

        # Sort SeqIDs: real sequences first (lexicographic), "." (globals) last
        all_seqids = sorted(
            all_seqids_seen.keys(),
            key=lambda s: (s == '.', s)
        )
        print(f'  Found {len(all_seqids)} unique SeqIDs across all samples.')

        # ── Phase 2: process each SeqID (parallel if threads > 1) ─────────────
        col_prefixes = [col_prefix for _, col_prefix in sample_info]
        worker_args = [
            (seqid,
             [ssp.get(seqid) for ssp in sample_seqid_paths],
             col_prefixes, temp_out_dir,
             metadata_cols, bp_meta, ALL_BPS, min_cov, decimals, report_non_qualified,
             min_samples_pct, min_group_pct)
            for seqid in all_seqids
        ]

        mode = f'{min(threads, len(all_seqids))} workers' if threads > 1 else 'sequential'
        print(f'Phase 2/3 — Processing {len(all_seqids)} SeqID chunks ({mode})...')

        if threads > 1:
            with multiprocessing.Pool(processes=min(threads, len(all_seqids))) as pool:
                results = pool.starmap(_process_seqid_chunk, worker_args)
        else:
            results = [_process_seqid_chunk(*a) for a in worker_args]

        # Restore stable SeqID order (parallel mode may return out of order)
        seqid_order = {s: i for i, s in enumerate(all_seqids)}
        results.sort(key=lambda r: seqid_order[r[0]])

        # ── Phase 3: concatenate per-SeqID chunks into final output files ─────────
        # Two output directories (one per metric type), each with 16 BP files.
        print('Phase 3/3 — Writing final output files...')
        output_files = []
        for metric in ('espf', 'espr'):
            out_dir = f'{output_prefix}_{metric}'
            os.makedirs(out_dir, exist_ok=True)
            for bp in ALL_BPS:
                bp_chunks = [
                    out_paths[(bp, metric)]
                    for _, out_paths in results
                    if (bp, metric) in out_paths
                ]
                if not bp_chunks:
                    continue

                out_path = os.path.join(out_dir, f'{bp}.tsv')
                with open(out_path, 'w') as fout:
                    with open(bp_chunks[0]) as first:
                        shutil.copyfileobj(first, fout)   # includes header
                    for chunk_path in bp_chunks[1:]:
                        with open(chunk_path) as f:
                            next(f)  # skip header
                            shutil.copyfileobj(f, fout)
                output_files.append(out_path)
                print(f'  {out_path}')

    print(f'\nDone. {len(sample_info)} samples, {len(all_seqids)} SeqIDs, '
          f'{len(output_files)} output files written '
          f'in {output_prefix}_espf/ and {output_prefix}_espr/.')


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
    min_samples_pct = None   # Default: no global-sample-% filter
    min_group_pct   = None   # Default: no per-group-% filter

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

        # Check for --min-samples-pct flag
        if arg.startswith('--min-samples-pct'):
            if '=' in arg:
                val = arg.split('=', 1)[1]
            else:
                try:
                    next_i = next(args_iter)
                    val = sys.argv[next_i]
                except StopIteration:
                    print('ERROR: --min-samples-pct requires a value', file=sys.stderr)
                    sys.exit(1)
            try:
                min_samples_pct = float(val)
            except ValueError:
                print('ERROR: --min-samples-pct requires a numeric value (0–100)', file=sys.stderr)
                sys.exit(1)
            if not (0.0 <= min_samples_pct <= 100.0):
                print('ERROR: --min-samples-pct must be between 0 and 100', file=sys.stderr)
                sys.exit(1)
            continue

        # Check for --min-group-pct flag
        if arg.startswith('--min-group-pct'):
            if '=' in arg:
                val = arg.split('=', 1)[1]
            else:
                try:
                    next_i = next(args_iter)
                    val = sys.argv[next_i]
                except StopIteration:
                    print('ERROR: --min-group-pct requires a value', file=sys.stderr)
                    sys.exit(1)
            try:
                min_group_pct = float(val)
            except ValueError:
                print('ERROR: --min-group-pct requires a numeric value (0–100)', file=sys.stderr)
                sys.exit(1)
            if not (0.0 <= min_group_pct <= 100.0):
                print('ERROR: --min-group-pct must be between 0 and 100', file=sys.stderr)
                sys.exit(1)
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
    result = merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id, min_cov, threads, decimals, report_non_qualified, min_samples_pct, min_group_pct)
    
    print("\nAnalysis complete!")