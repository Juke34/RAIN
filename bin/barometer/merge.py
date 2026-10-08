"""barometer.merge – Cross-mtype / cross-vtype "truly global" ranking.

Consumes the per-mtype ``global_ranking`` CSVs produced by the independent
``barometer_analyze`` runs (espf/aggregate, espf/feature, espf/sites,
espr/aggregate, espr/feature, espr/sites), merges them into a single
*truly global* ranking, and runs a second-pass analysis by redistributing
the significant biomarkers back to their value type.

The second pass is done per value type because ``analyze_section`` needs the
per-vtype sample columns (espf and espr have different column layouts), so a
cross-vtype second pass on a single merged DataFrame is not possible.

Raw inputs are provided as a per-vtype mapping (see
:func:`merge_global_rankings`), e.g.::

    {
      "espf": {"aggregates": "espf_agg.tsv", "features": "espf_feat.tsv",
               "sites": "espf_sites.tsv"},
      "espr": {"aggregates": "espr_agg.tsv", "features": "espr_feat.tsv",
               "sites": "espr_sites.tsv"}
    }
"""

import logging
import os

import pandas as pd

from .analysis import analyze_section
from .ranking import build_global_ranking
from .utils import (
    slurp_file,
    safe_mkdir,
    harmonize_columns,
    create_uid_column,
    parse_sample_columns,
    cols_for_vtype,
    sample_info_for_vtype,
    filter_significant,
)

log = logging.getLogger(__name__)

# Full second-pass test set (same as the per-vtype second pass in __main__).
SECOND_PASS_TESTS = [
    "descriptive", "multivariate", "correlation", "differential",
    "ranking", "classification", "stability", "heatmap",
]


def _load_raw_df(path):
    """Load a raw DRIP TSV and prepare it (strip, harmonize, uid, categories).

    Mirrors the preparation done in ``main()`` so that uids computed here match
    the uids present in the per-mtype global_ranking CSVs.
    """
    df = slurp_file(path, separator="\t")
    df.columns = df.columns.str.strip()
    df = harmonize_columns(df)
    df = create_uid_column(df)
    for col in df.select_dtypes(include=['object', 'string']).columns:
        df[col] = df[col].str.strip() if df[col].dtype in ['object', 'string'] else df[col]
    for col in ("Mtype", "Ptype", "Type", "Ctype", "Mode", "Strand", "SeqID"):
        if col in df.columns:
            df[col] = df[col].astype("category")
    return df


def _load_raw_pair(raw_inputs, vtype):
    """Load (agg_df, feat_df) for a main value type from the raw inputs mapping."""
    entry = raw_inputs.get(vtype, {})
    agg_path = entry.get("aggregates")
    feat_path = entry.get("features")
    agg_df = _load_raw_df(agg_path) if agg_path and os.path.exists(agg_path) else None
    feat_df = _load_raw_df(feat_path) if feat_path and os.path.exists(feat_path) else None
    return agg_df, feat_df


def merge_global_rankings(results_dirs, outdir, stat_test="auto",
                          raw_inputs=None, max_bmks=500, bmk_filter=None):
    """Merge the per-mtype global rankings into a truly-global ranking.

    Args:
        results_dirs: a single dir (str) or a list of dirs, each containing
            ``espf/`` and/or ``espr/`` subdirs holding the per-mtype
            (aggregate / feature / sites) ``global_ranking`` CSVs. In Nextflow
            each barometer_analyze invocation produces its own results dir, so
            several are passed and their CSVs are pooled.
        outdir: output dir for the merged ranking.
        stat_test: statistical test mode (beta-binomial / nonparametric / ...).
        raw_inputs: dict mapping vtype -> {"aggregates": path, "features": path,
            "sites": path}, used to rebuild the full rows for the second pass.
            May be None (second pass is then skipped).
        max_bmks: max BMKs for ML models.
        bmk_filter: cascade of *_padj columns for significant BMK selection.
    """
    safe_mkdir(outdir)
    raw_inputs = raw_inputs or {}
    if isinstance(results_dirs, str):
        results_dirs = [results_dirs]

    # ------------------------------------------------------------------
    # 1. Collect the per-mtype global_ranking_primary_padj.csv files
    #    Layout: <results_dir>/<vtype>/<mtype>/global_ranking/...
    # ------------------------------------------------------------------
    gr_csvs = []
    for results_dir in results_dirs:
        if not os.path.isdir(results_dir):
            continue
        for vtype in sorted(os.listdir(results_dir)):
            vtype_dir = os.path.join(results_dir, vtype)
            if not os.path.isdir(vtype_dir):
                continue
            for mtype in ("aggregate", "feature", "sites"):
                gr_csv = os.path.join(vtype_dir, mtype, "global_ranking",
                                      "global_ranking_primary_padj.csv")
                if os.path.isfile(gr_csv):
                    gr_csvs.append(gr_csv)

    if not gr_csvs:
        log.warning("No per-mtype global_ranking CSVs found in %s; nothing to merge", results_dirs)
        return outdir

    log.info(f"Merging {len(gr_csvs)} per-mtype global rankings...")
    merged = pd.concat([slurp_file(p, separator=",") for p in gr_csvs], ignore_index=True)
    log.info(f"  Merged table: {len(merged)} rows")

    # ------------------------------------------------------------------
    # 2. Build the truly-global ranking (reuses the per-mtype engine).
    #    build_global_ranking recomputes n_sections_total / n_sections_significant
    #    across ALL sections (all mtypes + vtypes) → the cross-validation metric
    #    now spans the whole dataset.
    # ------------------------------------------------------------------
    build_global_ranking(merged, outdir, stat_test)

    # ------------------------------------------------------------------
    # 3. Second pass by redistribution per value type
    # ------------------------------------------------------------------
    any_sig_csv = os.path.join(outdir, "significant_bmks_any_comparison", "significant_bmks.csv")
    if not os.path.isfile(any_sig_csv):
        log.warning("No merged significant BMKs; skipping second pass")
        return outdir

    df_sig_all = slurp_file(any_sig_csv, separator=",")
    log.info(f"  Second pass on {len(df_sig_all)} merged significant BMKs")

    for vtype in sorted(df_sig_all["value_type"].unique()):
        vtype_sig = df_sig_all[df_sig_all["value_type"] == vtype]

        # Route site BMKs (section == "sites_all") to the per-site raw data;
        # everything else uses the main per-vtype aggregates/features raw data.
        site_mask = vtype_sig["section"].astype(str) == "sites_all"
        main_sig = vtype_sig[~site_mask]
        site_sig = vtype_sig[site_mask]

        # --- Main BMKs (espf / espr) ---
        if len(main_sig) > 0:
            agg_df, feat_df = _load_raw_pair(raw_inputs, vtype)
            sample_df = agg_df if agg_df is not None else feat_df
            if sample_df is not None:
                sample_info = parse_sample_columns(sample_df.columns)
                vcols = cols_for_vtype(sample_info, vtype)
                v_sample_info = sample_info_for_vtype(sample_info, vtype)
                if vcols:
                    out = os.path.join(outdir, "second_pass", vtype)
                    safe_mkdir(out)
                    log.info(f"  Second pass {vtype} (main): {len(main_sig)} BMKs")
                    df_sig = filter_significant(main_sig, feat_df, agg_df, uid_col="uid",
                                                info=any_sig_csv)
                    if df_sig is not None and len(df_sig) > 0:
                        analyze_section(df_sig, vcols, v_sample_info, out, "Section Significant",
                                        stat_test=stat_test, bmk_filter_cols=bmk_filter,
                                        max_bmks=max_bmks, enabled_tests=SECOND_PASS_TESTS)

        # --- Site BMKs (per-site, per-vtype) ---
        if len(site_sig) > 0:
            sites_path = raw_inputs.get(vtype, {}).get("sites")
            sites_df = _load_raw_df(sites_path) if sites_path and os.path.exists(sites_path) else None
            if sites_df is not None:
                sites_sample_info = parse_sample_columns(sites_df.columns)
                vcols = cols_for_vtype(sites_sample_info, vtype)
                v_sample_info = sample_info_for_vtype(sites_sample_info, vtype)
                if vcols:
                    out = os.path.join(outdir, "second_pass", vtype, "sites")
                    safe_mkdir(out)
                    log.info(f"  Second pass {vtype} (sites): {len(site_sig)} BMKs")
                    df_sig = filter_significant(site_sig, sites_df, None, uid_col="uid",
                                                info=any_sig_csv)
                    if df_sig is not None and len(df_sig) > 0:
                        analyze_section(df_sig, vcols, v_sample_info, out, "Section Significant",
                                        stat_test=stat_test, bmk_filter_cols=bmk_filter,
                                        max_bmks=max_bmks, enabled_tests=SECOND_PASS_TESTS)

    return outdir
