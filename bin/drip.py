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

HELP_TEXT = """
DRIP - RNA Editing Analysis Tool

DESCRIPTION:
This script analyzes RNA editing from standardized pluviometer files. It calculates
two key metrics for all 12 genome-variant base pair combinations across multiple 
samples and combines them into a unified matrix format.
drip.py implicitly assumes that the combination of the 11 columns forms a unique key, used to merge rows across samples.

USAGE:
./drip.py --output OUTPUT_PREFIX FILE1:GROUP1:SAMPLE1:REP1 FILE2:GROUP2:SAMPLE2:REP2 [...]
./drip.py --help | -h

ARGUMENTS:
--output OUTPUT_PREFIX, -o OUTPUT_PREFIX
                 Prefix for the output directories (required).
                 Creates OUTPUT_PREFIX_espf/ and OUTPUT_PREFIX_espr/,
                 each containing 12 TSV files (AA.tsv, AC.tsv, …, TT.tsv).
FILEn:GROUPn:SAMPLEn:REPn  
                 Input file path, group name, sample name, and replicate ID
                 separated by colons. All four components are required.
--with-file-id   Include file ID in column names (default: omit file ID)
--preserve-covered-features
                 Include ::successes and ::trials count columns in the output
                 (for coverage/depth information).  Row filtering is unchanged:
                 the flag only adds count columns, it does not relax the
                 coverage/editing filters.  Intended for count-based differential
                 tests (e.g. beta-binomial) that need raw counts alongside
                 proportions.

Row filtering is done in two stages.  A row is kept only if it passes
BOTH stages:

    keep = COVERAGE AND EDITING

    COVERAGE = (min-samples AND min-samples-pct)
              OR (min-group-samples AND min-group-samples-pct)

    EDITING  = (min-group-samples-edited AND min-group-samples-pct-edited)

Definitions:
  - "covered" = non-NA (0.0 counts as covered: the position was observed,
    just without editing).
  - "edited"  = non-NA AND non-zero (an actual editing event).

COVERAGE filters (count covered, 0.0 included):
  Flag                        Counter       AND with                  OR with
  --min-samples N             all samples   --min-samples-pct         group criterion
  --min-samples-pct X         all samples   --min-samples             group criterion
  --min-group-samples N       per group     --min-group-samples-pct   global criterion
  --min-group-samples-pct Y   per group     --min-group-samples       global criterion

  - "global criterion" = (min-samples AND min-samples-pct) evaluated on
    the total number of samples across all groups.
  - "group criterion"  = (min-group-samples AND min-group-samples-pct)
    evaluated per group; the row passes if ANY single group satisfies
    both conditions of the pair.
  - If only one flag of an AND pair is given, only that flag applies.
  - If no coverage flag is given, the coverage stage is not applied
    (all rows pass).
  - --min-group-samples N guards against single-sample noise in small
    groups (where 1/1 = 100% would pass --min-group-samples-pct alone).

EDITING filters (count edited, non-NA and non-zero):
  Flag                                Counter    AND with                        OR with
  --min-group-samples-edited N        per group  --min-group-samples-pct-edited  other groups
  --min-group-samples-pct-edited Y    per group  --min-group-samples-edited      other groups

  - Both conditions of the pair must be met by the SAME group; the row
    passes if ANY single group satisfies them.
  - Defaults: both None (disabled) → the editing stage is not applied.
    Set N=1 to require at least one edited sample in some group.

Rules:
  - 0.0 counts as covered but NOT as edited.
  - With no filter flag set, all rows are kept (including rows where
    every cell is NA).
--min-cov N      Minimum read coverage threshold (default: 1). Positions with a
                 denominator (genome base count for espf, read count for espr)
                 strictly below this value are reported as NA instead of 0.
--threads N, -t N
                 Number of parallel threads to use for writing output files
                 (default: 1, sequential). Max useful value is 12 (one per base pair).
--decimals N, -d N
                 Number of decimal places for output values (default: 4).
                 Reduces file size by rounding espf/espr metrics.
--bps LIST       Comma-separated base pairs to analyze, e.g. 'AC,AG,CA'.
                 Only the selected BPs are computed and written (one output
                 file per BP per metric). Default: all 12
                 (AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG).
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
For each line, the script calculates metrics for all 12 base pair combinations:

For each combination XY (where X = genome base, Y = read base):

1. XY_espf (edited_sites_proportion_feature) - Proportion of XY sites in the DNA feature:
   Formula: XY_SiteBasePairingsQualified / X_QualifiedBases
   This represents the proportion of qualified X positions that show X-to-Y variation in the feature.
     Count columns preserve XY_SiteBasePairingsQualified as successes and X_QualifiedBases as trials.
     In other terms for A-to-G example:
     - successes = number of reference sites A where a G edit was observed;
     - trials = total number of reference sites A qualified in the feature, whether they harbor editing or not.

2. XY_espr (edited_sites_proportion_reads) - Proportion of XY pairing in reads:
   Formula: XY_ReadBasePairings / (XA + XC + XG + XT)_ReadBasePairings
   This represents the proportion of X-position reads that show Y in the reads.
     Count columns preserve XY_ReadBasePairings as successes and the X-read sum as trials.
     In other terms for A-to-G example:
        - successes = reads that cover a position where the reference is A and that observe G;
        - trials = all reads covering these positions where the reference is A, regardless of whether they observe A, C, G, or T.
     
All 12 combinations are calculated: AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG

OUTPUT FORMAT:
Two output directories, one per metric type, each containing 12 TSV files
(one per base pair combination).  Row filtering (--min-samples-pct, etc.)
is applied independently per metric: a row may appear in the espf output
but not in the espr output and vice versa.

Directories created:
- OUTPUT_PREFIX_espf/   espf metric (proportion of edited sites in feature)
- OUTPUT_PREFIX_espr/   espr metric (proportion of edited reads)

Files in each directory (12 per metric):
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
- GROUP::SAMPLE::REPLICATE::<metric>::successes and ::trials for count-based models

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

def print_help():
    """Print help message explaining the script usage and calculations."""
    print(HELP_TEXT)
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
    gb_offset = gb_idx * 4               # first XA/XC/XG/XT index in 12-value vector

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
    result[f'{col_prefix}::espf::successes'] = bp_sites.where(mask_f)
    result[f'{col_prefix}::espf::trials'] = x_count.where(mask_f)
    result[f'{col_prefix}::espf'] = np.where(
        mask_f, bp_sites / x_count.where(mask_f, 1), np.nan
    )
    result[f'{col_prefix}::espf'] = result[f'{col_prefix}::espf'].round(decimals)

    mask_r = total_reads >= min_cov
    result[f'{col_prefix}::espr::successes'] = bp_reads.where(mask_r)
    result[f'{col_prefix}::espr::trials'] = total_reads.where(mask_r)
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
    result[f'{col_prefix}::espf::successes'] = bp_sites.where(mask_f)
    result[f'{col_prefix}::espf::trials'] = x_count.where(mask_f)
    result[f'{col_prefix}::espf'] = np.where(
        mask_f, bp_sites / x_count.where(mask_f, 1), np.nan
    ).round(decimals)

    mask_r = total_reads >= min_cov
    result[f'{col_prefix}::espr::successes'] = bp_reads.where(mask_r)
    result[f'{col_prefix}::espr::trials'] = total_reads.where(mask_r)
    result[f'{col_prefix}::espr'] = np.where(
        mask_r, bp_reads / total_reads.where(mask_r, 1), np.nan
    ).round(decimals)

    return result


def _write_one_bp(bp, accumulated, output_prefix, metadata_cols, report_non_qualified):
    """Sort, filter and write the accumulated DataFrame for one BP."""
    merged = accumulated.sort_values(['SeqID', 'ParentIDs', 'Mode'])

    if not report_non_qualified:
        metric_cols = [c for c in merged.columns
                       if c not in metadata_cols and c.endswith(('::espf', '::espr'))]
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
                         min_samples=None, min_samples_pct=None,
                         min_group_samples_pct=None, min_group_samples=None,
                         min_group_samples_edited=None, min_group_samples_pct_edited=None,
                         preserve_covered_features=False):
    """Compute the selected BPs for one SeqID and write per-BP chunk files.

    Each element in chunk_paths_by_sample is a path (str) or None when that
    sample has no rows for this SeqID — the outer join fills NaN for it.
    Returns (seqid, {bp: out_path, ...}) — only BPs with at least one row.
    """
    mixed_dtypes_local = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}
    bp_accumulators = [None] * len(all_bps)

    for col_prefix, chunk_path in zip(col_prefixes, chunk_paths_by_sample):
        if chunk_path is None:
            # Crée un DataFrame vide avec les bonnes colonnes de métadonnées et une colonne NA pour ce sample
            # On suppose qu'il y a au moins un sample avec données pour récupérer la structure
            # On construit une ligne unique avec NA pour les colonnes metrics, et les métadonnées à la valeur du seqid
            # Pour chaque BP, on crée une ligne NA
            for i in range(len(all_bps)):
                # On crée un DataFrame avec une seule ligne, toutes les métadonnées à seqid ou '.' et la colonne metric à NA
                meta = {k: seqid if k == 'SeqID' else '.' for k in metadata_cols}
                # Une seule ligne
                row = {
                    **meta,
                    f'{col_prefix}::espf': np.nan,
                    f'{col_prefix}::espf::successes': np.nan,
                    f'{col_prefix}::espf::trials': np.nan,
                    f'{col_prefix}::espr': np.nan,
                    f'{col_prefix}::espr::successes': np.nan,
                    f'{col_prefix}::espr::trials': np.nan,
                }
                bp_data = pd.DataFrame([row])
                if bp_accumulators[i] is None:
                    bp_accumulators[i] = bp_data
                else:
                    bp_accumulators[i] = bp_accumulators[i].merge(
                        bp_data, on=metadata_cols, how='outer'
                    )
                    del bp_data
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
            # --preserve-covered-features only controls whether the raw count
            # columns (::successes / ::trials) are written to the output, so
            # downstream count-based tests (e.g. beta-binomial) can see the
            # observation depth behind each proportion (including 0.0 values).
            # It does NOT affect row filtering.
            if preserve_covered_features:
                count_cols = [c for c in acc.columns
                              if c not in metadata_cols
                              and c.endswith((f'::{metric}::successes', f'::{metric}::trials'))]
            else:
                count_cols = []
            acc_metric = acc[metadata_cols + metric_cols + count_cols]

            # Cell-level masks:
            #   covered = non-NA.  NA means the position was not covered
            #     (or below min_cov, or absent from this sample via the
            #     outer join).  0.0 means covered but no editing observed.
            #   edited  = non-NA AND non-zero (an actual editing event).
            covered = acc_metric[metric_cols].notna()
            edited = covered & (acc_metric[metric_cols] != 0)

            # Row-level decision, applied independently per metric type
            # (espf rows and espr rows are filtered separately):
            #
            #   keep = COVERAGE AND EDITING
            #
            #   COVERAGE = (min-samples AND min-samples-pct)
            #               OR (min-group-samples AND min-group-samples-pct)
            #     - global criterion: both flags AND on the total covered
            #       sample count (all groups).
            #     - group criterion: both flags AND on the covered count
            #       of ONE group; OR across groups.  Groups are identified
            #       by the first :: component of the column name.
            #     - one flag of a pair given → only that flag applies.
            #     - no coverage flag → coverage stage not applied (all
            #       rows pass).
            #
            #   EDITING = (min-group-samples-edited
            #               AND min-group-samples-pct-edited)
            #     - both flags AND on the edited count of ONE group;
            #       OR across groups.
            #     - defaults (None, None) → editing stage not applied.
            groups: dict[str, list[str]] = {}
            for c in metric_cols:
                groups.setdefault(c.split('::')[0], []).append(c)
            if (min_samples is None and min_samples_pct is None
                    and min_group_samples_pct is None and min_group_samples is None):
                keep_cov = pd.Series(True, index=acc_metric.index)
            else:
                keep_cov = pd.Series(False, index=acc_metric.index)
                if min_samples is not None or min_samples_pct is not None:
                    n = len(metric_cols)
                    n_cov = covered.sum(axis=1)
                    samples_ok = pd.Series(True, index=acc_metric.index)
                    if min_samples is not None:
                        samples_ok &= n_cov >= min_samples
                    if min_samples_pct is not None:
                        samples_ok &= n_cov / n >= min_samples_pct / 100.0
                    keep_cov |= samples_ok
                if (min_group_samples_pct is not None
                        or min_group_samples is not None):
                    for g_cols in groups.values():
                        n_cov = covered[g_cols].sum(axis=1)
                        group_ok = pd.Series(True, index=acc_metric.index)
                        if min_group_samples is not None:
                            group_ok &= n_cov >= min_group_samples
                        if min_group_samples_pct is not None:
                            group_ok &= (
                                n_cov / len(g_cols) >= min_group_samples_pct / 100.0
                            )
                        keep_cov |= group_ok

            # Stage 2 — editing.  Skipped entirely when both flags are None
            # (the default), so the default behavior is coverage-only.
            if (min_group_samples_edited is None
                    and min_group_samples_pct_edited is None):
                keep = keep_cov
            else:
                keep_edited = pd.Series(False, index=acc_metric.index)
                for g_cols in groups.values():
                    n_ed = edited[g_cols].sum(axis=1)
                    group_edited_ok = pd.Series(True, index=acc_metric.index)
                    if min_group_samples_edited is not None:
                        group_edited_ok &= n_ed >= min_group_samples_edited
                    if min_group_samples_pct_edited is not None:
                        group_edited_ok &= (
                            n_ed / len(g_cols) >= min_group_samples_pct_edited / 100.0
                        )
                    keep_edited |= group_edited_ok
                keep = keep_cov & keep_edited

            acc_metric = acc_metric[keep]

            if len(acc_metric) > 0:
                out_path = os.path.join(temp_out_dir, metric, f'{metric}_{bp}_{safe}.tsv')
                acc_metric.to_csv(out_path, sep='\t', index=False, na_rep='NA')
                out_paths[(bp, metric)] = out_path

        bp_accumulators[i] = None

    gc.collect()
    return seqid, out_paths


def _output_header(metadata_cols, col_prefixes, metric, preserve_covered_features):
    """Column names of a per-BP output file.

    Mirrors the order built in _process_seqid_chunk: metadata, then one
    metric column per sample, then (when preserve_covered_features) the
    ::successes and ::trials columns per sample.
    """
    cols = list(metadata_cols)
    for col_prefix in col_prefixes:
        cols.append(f'{col_prefix}::{metric}')
    if preserve_covered_features:
        for col_prefix in col_prefixes:
            cols.append(f'{col_prefix}::{metric}::successes')
            cols.append(f'{col_prefix}::{metric}::trials')
    return cols


def merge_samples(file_group_sample_replicate_dict, output_prefix, include_file_id=False,
                  min_cov=1, threads=1, decimals=4,
                  min_samples=None, min_samples_pct=None, min_group_samples_pct=None,
                  min_group_samples=None, min_group_samples_edited=None,
                  min_group_samples_pct_edited=None, preserve_covered_features=False,
                  bps=None):
    """Produce one output file per selected base pair combination.

    Memory strategy — three phases:
      Phase 1 (Split):   each input file is read once and split by SeqID into
                         temporary chunk files.  Peak RAM = one full input file.
      Phase 2 (Process): each SeqID is processed independently (all 12 BPs,
                         across all samples) and results written to temp chunks.
                         Peak RAM per worker ≈ num_samples × one-SeqID slice.
                         Workers run in parallel when threads > 1.
      Phase 3 (Concat):  per-SeqID temp chunks are appended in order to the
                         12 final output files.  Peak RAM = one chunk at a time.
    """
    ALL_BPS = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT',
               'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
    ALL_PAIRINGS = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                    'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    BASES = ['A', 'C', 'G', 'T']
    if bps is None:
        bps = ALL_BPS
    else:
        invalid = [bp for bp in bps if bp not in ALL_BPS]
        if invalid:
            raise ValueError(f"Invalid base pair(s): {', '.join(invalid)}. "
                             f"Valid values: {', '.join(ALL_BPS)}")
        if not bps:
            raise ValueError("At least one base pair must be selected")
    metadata_cols = ['SeqID', 'ParentIDs', 'ID', 'Mtype', 'Ptype', 'Type',
                     'Ctype', 'Mode', 'Start', 'End', 'Strand']
    needed_cols = metadata_cols + ['QualifiedBases',
                                   'SiteBasePairingsQualified',
                                   'ReadBasePairingsQualified']
    mixed_dtypes = {'SeqID': str, 'Start': str, 'End': str, 'Strand': str}
    bp_meta = [(ALL_PAIRINGS.index(bp), BASES.index(bp[0]), BASES.index(bp[0]) * 4)
               for bp in bps]

    sample_info = []
    for filepath, (group_name, sample_name, replicate) in file_group_sample_replicate_dict.items():
        filename_stem = Path(filepath).stem
        file_id = '_'.join(filename_stem.split('_')[:-1])
        col_prefix = (f'{group_name}::{sample_name}::{replicate}::{file_id}'
                      if include_file_id
                      else f'{group_name}::{sample_name}::{replicate}')
        sample_info.append((filepath, col_prefix))

    col_prefixes = [col_prefix for _, col_prefix in sample_info]

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

        if len(all_seqids) == 0:
            print('  WARNING: No data rows found across input files (headers only or empty content).')
            # Write header-only files so downstream consumers (Nextflow)
            # still find the expected output files.
            for metric in ('espf', 'espr'):
                out_dir = f'{output_prefix}_{metric}'
                os.makedirs(out_dir, exist_ok=True)
                for bp in bps:
                    out_path = os.path.join(out_dir, f'{output_prefix}_{metric}_{bp}.tsv')
                    with open(out_path, 'w') as fout:
                        fout.write('\t'.join(_output_header(
                            metadata_cols, col_prefixes, metric,
                            preserve_covered_features)) + '\n')
            print(f'\nDone. {len(sample_info)} samples, 0 SeqIDs, 0 data rows '
                  f'(header-only files written) in {output_prefix}_espf/ and '
                  f'{output_prefix}_espr/.')
            return []

        # ── Phase 2: process each SeqID (parallel if threads > 1) ─────────────
        worker_args = [
            (seqid,
             [ssp.get(seqid) for ssp in sample_seqid_paths],
             col_prefixes, temp_out_dir,
             metadata_cols, bp_meta, bps, min_cov, decimals,
             min_samples, min_samples_pct, min_group_samples_pct, min_group_samples,
             min_group_samples_edited, min_group_samples_pct_edited,
             preserve_covered_features)
            for seqid in all_seqids
        ]

        mode = f'{min(threads, len(all_seqids))} workers' if threads > 1 else 'sequential'
        print(f'Phase 2/3 — Processing {len(all_seqids)} SeqID chunks ({mode})...')

        if threads > 1 and len(all_seqids) > 0:
            with multiprocessing.Pool(processes=min(threads, len(all_seqids))) as pool:
                results = pool.starmap(_process_seqid_chunk, worker_args)
        else:
            results = [_process_seqid_chunk(*a) for a in worker_args]

        # Restore stable SeqID order (parallel mode may return out of order)
        seqid_order = {s: i for i, s in enumerate(all_seqids)}
        results.sort(key=lambda r: seqid_order[r[0]])

        # ── Phase 3: concatenate per-SeqID chunks into final output files ─────────
        # Two output directories (one per metric type), each with 12 BP files.
        print('Phase 3/3 — Writing final output files...')
        output_files = []
        for metric in ('espf', 'espr'):
            out_dir = f'{output_prefix}_{metric}'
            os.makedirs(out_dir, exist_ok=True)
            for bp in bps:
                bp_chunks = [
                    out_paths[(bp, metric)]
                    for _, out_paths in results
                    if (bp, metric) in out_paths
                ]
                out_path = os.path.join(out_dir, f'{output_prefix}_{metric}_{bp}.tsv')
                if not bp_chunks:
                    # No rows survived the filters for this BP: write a
                    # header-only file so downstream consumers (Nextflow)
                    # still find the expected output file.
                    with open(out_path, 'w') as fout:
                        fout.write('\t'.join(_output_header(
                            metadata_cols, col_prefixes, metric,
                            preserve_covered_features)) + '\n')
                else:
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
    import argparse

    # Check for help flag
    if len(sys.argv) > 1 and sys.argv[1] in ['--help', '-h', 'help']:
        print_help()

    def _nonneg_int(value):
        ivalue = int(value)
        if ivalue < 0:
            raise argparse.ArgumentTypeError(f"must be >= 0, got {ivalue}")
        return ivalue

    def _pct(value):
        fvalue = float(value)
        if not 0.0 <= fvalue <= 100.0:
            raise argparse.ArgumentTypeError(f"must be between 0 and 100, got {fvalue}")
        return fvalue

    def _bps(value):
        valid = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT',
                 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
        parsed = [bp.strip().upper() for bp in value.split(',') if bp.strip()]
        if not parsed:
            raise argparse.ArgumentTypeError("at least one base pair is required")
        invalid = [bp for bp in parsed if bp not in valid]
        if invalid:
            raise argparse.ArgumentTypeError(
                f"invalid base pair(s): {', '.join(invalid)}. "
                f"Valid values: {', '.join(valid)}")
        return parsed

    parser = argparse.ArgumentParser(
        prog="drip.py",
        description=HELP_TEXT,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "samples", nargs="+", metavar="FILE:GROUP:SAMPLE:REPLICATE",
        help="One or more sample files. Each entry is "
             "FILE:GROUP:SAMPLE:REPLICATE (4 colon-separated parts). "
             "FILE must exist. Duplicate GROUP:SAMPLE:REPLICATE keys are "
             "renamed with a numeric suffix on the replicate.",
    )
    parser.add_argument("-o", "--output", required=True, metavar="PREFIX",
                        help="Output prefix (required).")
    parser.add_argument("--with-file-id", action="store_true",
                        help="Include the file name as an extra column in the output.")
    parser.add_argument("--preserve-covered-features", action="store_true",
                        help="Keep covered features and add coverage columns.")
    parser.add_argument("--min-cov", type=int, default=1, metavar="N",
                        help="Minimum coverage (default: 1).")
    parser.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                        help="Number of threads (default: 1).")
    parser.add_argument("-d", "--decimals", type=int, default=4, metavar="N",
                        help="Number of decimals in the output (default: 4).")
    parser.add_argument("--min-samples", type=_nonneg_int, default=None, metavar="N",
                        help="Coverage: minimum number of covered samples (AND with --min-samples-pct).")
    parser.add_argument("--min-samples-pct", type=_pct, default=None, metavar="PCT",
                        help="Coverage: minimum %% of covered samples (AND with --min-samples).")
    parser.add_argument("--min-group-samples", type=_nonneg_int, default=None, metavar="N",
                        help="Coverage: minimum covered samples in one group (AND with --min-group-samples-pct).")
    parser.add_argument("--min-group-samples-pct", type=_pct, default=None, metavar="PCT",
                        help="Coverage: minimum %% of covered samples in one group (AND with --min-group-samples).")
    parser.add_argument("--min-group-samples-edited", type=_nonneg_int, default=None, metavar="N",
                        help="Editing: minimum edited samples in one group (AND with --min-group-samples-pct-edited). Default: None (editing stage disabled).")
    parser.add_argument("--min-group-samples-pct-edited", type=_pct, default=None, metavar="PCT",
                        help="Editing: minimum %% of edited samples in one group (AND with --min-group-samples-edited). Default: None (editing stage disabled).")
    parser.add_argument("--bps", type=_bps, default=None, metavar="LIST",
                        help="Comma-separated base pairs to analyze, e.g. 'AC,AG,CA'. "
                             "Only the selected BPs are computed and written. "
                             "Default: all 12 (AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG).")
    args = parser.parse_args()

    # Validate threads
    if args.threads < 1:
        print("Error: --threads must be at least 1")
        sys.exit(1)
    if args.decimals < 0:
        print("Error: --decimals must be non-negative")
        sys.exit(1)

    # Parse sample files: FILE:GROUP:SAMPLE:REPLICATE
    file_group_sample_replicate_dict = {}
    group_sample_rep_counts = {}
    for sample_file in args.samples:
        parts = sample_file.split(':')
        if len(parts) != 4:
            print(f"ERROR: Invalid argument format '{sample_file}'", file=sys.stderr)
            print("Expected format: FILE:GROUP:SAMPLE:REPLICATE (all 4 components required)", file=sys.stderr)
            sys.exit(1)
        file_path, group, sample, replicate = parts
        if not os.path.exists(file_path):
            print(f"ERROR: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)
        # Handle duplicate group:sample:replicate combinations by adding a suffix
        group_sample_rep_key = f"{group}:{sample}:{replicate}"
        if group_sample_rep_key in group_sample_rep_counts:
            group_sample_rep_counts[group_sample_rep_key] += 1
            replicate = f"{replicate}_{group_sample_rep_counts[group_sample_rep_key]}"
            print(f"WARNING: Duplicate group:sample:replicate '{group_sample_rep_key}' found. "
                  f"Renaming replicate to '{replicate}'", file=sys.stderr)
        else:
            group_sample_rep_counts[group_sample_rep_key] = 1
        file_group_sample_replicate_dict[file_path] = (group, sample, replicate)

    output_prefix = args.output
    include_file_id = args.with_file_id
    min_cov = args.min_cov
    threads = args.threads
    decimals = args.decimals
    min_samples = args.min_samples
    min_samples_pct = args.min_samples_pct
    min_group_samples_pct = args.min_group_samples_pct
    min_group_samples = args.min_group_samples
    min_group_samples_edited = args.min_group_samples_edited
    min_group_samples_pct_edited = args.min_group_samples_pct_edited
    preserve_covered_features = args.preserve_covered_features


    # Process all samples
    result = merge_samples(
        file_group_sample_replicate_dict, output_prefix, include_file_id, min_cov,
        threads, decimals, min_samples, min_samples_pct,
        min_group_samples_pct, min_group_samples, min_group_samples_edited,
        min_group_samples_pct_edited, preserve_covered_features,
        bps=args.bps
    )
    
    print("\nAnalysis complete!")