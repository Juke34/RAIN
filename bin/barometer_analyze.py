#!/usr/bin/env python3
"""
barometer_analyze.py – Exhaustive analysis pipeline for rain biomarker data.

Loads barometer_aggregates_AG.tsv and barometer_features_AG.tsv, performs QC,
descriptive statistics, differential editing analysis, multivariate analysis,
correlation/network analysis, feature selection / biomarker ranking,
classification, stability analysis, and generates all tables + matplotlib figures.

Results are saved into an output directory (default: barometer_results/) structured
by value_type → mtype → section.

Usage:
    python barometer_analyze.py [--aggregates FILE] [--features FILE] [--outdir DIR]
"""

import argparse
import gc
import json
import logging
import multiprocessing
import os
import sys
import time
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

try:
    import psutil
except ImportError:
    psutil = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_mkdir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def prepare_df_for_task(df):
    """Convert DataFrame to a memory-efficient format for task submission.
    
    Uses pickle which is much faster than .to_dict('list') for large DataFrames.
    Returns a tuple ('pickled', bytes_data) that can be unpickled in the worker.
    """
    import pickle
    return ('pickled', pickle.dumps(df, protocol=pickle.HIGHEST_PROTOCOL))


def parse_sample_columns(columns):
    """Parse sample column names like 'test1::rain_chr21_small::rep1::espf'.

    Returns a list of dicts with keys: col, group, sample, rep, value_type.
    """
    parsed = []
    meta_cols = {"SeqID", "ParentIDs", "ID", "Mtype", "Ptype", "Type", "Ctype", "Mode", "Start", "End", "Strand"}
    for c in columns:
        if c in meta_cols:
            continue
        parts = c.split("::")
        if len(parts) == 4:
            parsed.append({
                "col": c,
                "group": parts[0],
                "sample": parts[1],
                "rep": parts[2],
                "value_type": parts[3],
            })
    return parsed


def get_value_types(sample_info):
    return sorted(set(s["value_type"] for s in sample_info))


def cols_for_vtype(sample_info, vtype):
    return [s["col"] for s in sample_info if s["value_type"] == vtype]


def group_for_col(sample_info, col):
    for s in sample_info:
        if s["col"] == col:
            return s["group"]
    return None


def sample_info_for_vtype(sample_info, vtype):
    return [s for s in sample_info if s["value_type"] == vtype]


def numeric_df(df, cols):
    """Return a copy with the given columns cast to numeric (coerce errors)."""
    out = df.copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def save_fig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Quality Control
# ---------------------------------------------------------------------------

def qc_analysis(df, sample_cols, outdir):
    """Basic quality control: missing values, distributions, outliers."""
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # Missing value counts
    log.info("Missing value counts")
    missing = ndf[sample_cols].isnull().sum()
    total = len(ndf)
    missing_pct = (missing / total * 100).round(2)
    miss_df = pd.DataFrame({"missing_count": missing, "missing_pct": missing_pct})
    miss_df.to_csv(os.path.join(outdir, "missing_values.csv"))
    results["missing"] = miss_df.to_dict()

    # Distribution summary per sample
    log.info("Distribution summary per sample")
    desc = ndf[sample_cols].describe().T
    desc.to_csv(os.path.join(outdir, "distribution_summary.csv"))

    # Box plot of all samples
    log.info("Box plot of all samples")
    fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.8), 5))
    ndf[sample_cols].boxplot(ax=ax, rot=90)
    ax.set_title("Distribution per sample")
    ax.set_ylabel("Value")
    save_fig(fig, os.path.join(outdir, "boxplot_samples.png"))

    # Heatmap of missing values (limit to avoid matplotlib memory errors with large datasets)
    log.info("Heatmap of missing values")
    max_heatmap_rows = 10000  # limit to prevent memory issues
    heatmap_data = ndf[sample_cols].isnull().astype(int)
    
    if len(heatmap_data) > max_heatmap_rows:
        # Sample: prioritize BMKs with most missing values
        missing_counts = heatmap_data.sum(axis=1)
        top_missing_idx = missing_counts.nlargest(max_heatmap_rows).index
        heatmap_data = heatmap_data.loc[top_missing_idx]
        log.info(f"  Limiting missing values heatmap to {max_heatmap_rows} BMKs with most missing (from {len(ndf)} total)")
    
    fig_height = min(20, max(4, len(heatmap_data) * 0.02))  # Cap at 20 inches
    try:
        fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.6), fig_height))
        sns.heatmap(heatmap_data, cbar=False, ax=ax, yticklabels=False)
        ax.set_title(f"Missing values heatmap ({len(heatmap_data)} BMKs)")
        save_fig(fig, os.path.join(outdir, "missing_heatmap.png"))
        log.info("  Missing values heatmap saved")
    except Exception as e:
        log.warning(f"  Missing values heatmap failed: {e}")

    results["desc"] = desc.to_dict()
    return results


# ---------------------------------------------------------------------------
# Descriptive Statistics
# ---------------------------------------------------------------------------

def descriptive_stats(df, sample_cols, sample_info, outdir):
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # Per-group statistics
    groups = sorted(set(s["group"] for s in sample_info))
    group_stats = {}
    for g in groups:
        gcols = [s["col"] for s in sample_info if s["group"] == g]
        vals = ndf[gcols].values.flatten()
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0:
            continue
        group_stats[g] = {
            "n": int(len(vals)),
            "mean": float(np.mean(vals)),
            "median": float(np.median(vals)),
            "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
            "q25": float(np.percentile(vals, 25)),
            "q75": float(np.percentile(vals, 75)),
        }
    gs_df = pd.DataFrame(group_stats).T
    gs_df.to_csv(os.path.join(outdir, "group_statistics.csv"))
    results["group_stats"] = group_stats

    # Per-BMK mean by group
    bmk_means = {}
    for g in groups:
        gcols = [s["col"] for s in sample_info if s["group"] == g]
        bmk_means[g] = ndf[gcols].mean(axis=1)
    bmk_mean_df = pd.DataFrame(bmk_means, index=df.index)
    if "ID" in df.columns:
        bmk_mean_df.index = df["ID"]
    bmk_mean_df.to_csv(os.path.join(outdir, "bmk_mean_by_group.csv"))

    # Violin plot by group (optimized with pd.melt)
    sample_cols_list = [s["col"] for s in sample_info]
    if sample_cols_list:
        # Create group mapping for columns
        col_to_group = {s["col"]: s["group"] for s in sample_info}
        # Melt the dataframe efficiently
        melted = ndf[sample_cols_list].melt(var_name="sample_col", value_name="value")
        melted["group"] = melted["sample_col"].map(col_to_group)
        melted = melted.dropna(subset=["value"])
        
        if len(melted) > 0:
            fig, ax = plt.subplots(figsize=(max(6, len(groups) * 2), 5))
            sns.violinplot(data=melted, x="group", y="value", ax=ax, inner="box")
            ax.set_title("Value distribution by group")
            save_fig(fig, os.path.join(outdir, "violin_by_group.png"))

    return results


# ---------------------------------------------------------------------------
# Differential Editing Analysis
# ---------------------------------------------------------------------------

def test_normality_and_homogeneity(group_values):
    """Test for normality (Shapiro-Wilk) and homoscedasticity (Bartlett).
    
    Returns:
        dict with keys: is_normal, is_homogeneous, shapiro_pvals, bartlett_pval
    """
    results = {
        "is_normal": False,
        "is_homogeneous": False,
        "shapiro_pvals": [],
        "bartlett_pval": None
    }
    
    non_empty = [gv for gv in group_values.values() if len(gv) >= 3]  # Shapiro needs n>=3
    if len(non_empty) < 2:
        return results
    
    # Test normality per group (Shapiro-Wilk)
    shapiro_pvals = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        for gv in non_empty:
            if len(gv) >= 3:
                try:
                    # Check for near-constant data before testing
                    if np.std(gv) < 1e-10:
                        continue  # Skip constant data
                    _, p = stats.shapiro(gv)
                    shapiro_pvals.append(p)
                except Exception:
                    pass
    
    results["shapiro_pvals"] = shapiro_pvals
    results["is_normal"] = all(p > 0.05 for p in shapiro_pvals) if shapiro_pvals else False
    
    # Test homogeneity of variances (Bartlett)
    if len(non_empty) >= 2:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            try:
                # Check all groups have variance
                if all(np.std(gv) > 1e-10 for gv in non_empty):
                    _, p = stats.bartlett(*non_empty)
                    results["bartlett_pval"] = p
                    results["is_homogeneous"] = p > 0.05
            except Exception:
                pass
    
    return results


def differential_analysis(df, sample_cols, sample_info, outdir, stat_test="auto"):
    """Pairwise and global differential tests between groups.
    
    Args:
        stat_test: 'auto', 'parametric', 'nonparametric', 'welch', 'kruskal'
            - auto: test normality/homogeneity and choose best test
            - parametric: Student t-test / ANOVA (assumes normality + equal variances)
            - nonparametric: Mann-Whitney U / Kruskal-Wallis (no assumptions)
            - welch: Welch t-test / Welch ANOVA (assumes normality, unequal variances OK)
            - kruskal: alias for nonparametric (backward compatibility)
    
    Note: Expects df to already be pre-filtered for variable BMKs (done in analyze_section)
    """
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    groups = sorted(set(s["group"] for s in sample_info))
    results = {"pairwise": {}, "global": {}, "stat_test_used": stat_test}

    if len(groups) < 2:
        log.warning("Less than 2 groups, skipping differential analysis.")
        return results
    
    if len(ndf) == 0:
        log.warning("No BMKs remaining for differential analysis")
        return results
    
    # Normalize test names
    if stat_test == "kruskal":
        stat_test = "nonparametric"

    # OPTIMIZATION: Pre-compute group columns to avoid O(N_samples) per BMK per group
    group_cols = {g: [s["col"] for s in sample_info if s["group"] == g] for g in groups}

    rows = []
    
    # Suppress scipy warnings about numerical issues (we handle them with try/except)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Precision loss occurred")
        warnings.filterwarnings("ignore", message="invalid value encountered")
        warnings.filterwarnings("ignore", message="divide by zero encountered")
        
        for idx in ndf.index:
            row_data = {"index": idx}
            if "ID" in df.columns:
                row_data["ID"] = df.loc[idx, "ID"]

            group_values = {}
            for g in groups:
                gcols = group_cols[g]  # Use pre-computed dict
                vals = ndf.loc[idx, gcols].dropna().values.astype(float)
                group_values[g] = vals
                row_data[f"mean_{g}"] = np.mean(vals) if len(vals) > 0 else np.nan
                row_data[f"n_{g}"] = len(vals)

            non_empty = [gv for gv in group_values.values() if len(gv) > 0]
            
            # Determine which test to use
            test_choice = stat_test
            if stat_test == "auto" and len(non_empty) >= 2:
                # Test assumptions
                assumptions = test_normality_and_homogeneity(group_values)
                row_data["shapiro_pvals"] = ",".join([f"{p:.4f}" for p in assumptions["shapiro_pvals"]]) if assumptions["shapiro_pvals"] else ""
                row_data["bartlett_pval"] = assumptions["bartlett_pval"]
                
                # Choose test based on assumptions
                if assumptions["is_normal"] and assumptions["is_homogeneous"]:
                    test_choice = "parametric"
                elif assumptions["is_normal"] and not assumptions["is_homogeneous"]:
                    test_choice = "welch"
                else:
                    test_choice = "nonparametric"
                
                row_data["test_selected"] = test_choice
            
            # Global tests (comparing all groups)
            if len(non_empty) >= 2 and all(len(v) >= 1 for v in non_empty):
                # Kruskal-Wallis (always compute for backward compatibility)
                try:
                    stat, pval = stats.kruskal(*non_empty)
                    row_data["kruskal_stat"] = stat
                    row_data["kruskal_pval"] = pval
                except Exception:
                    row_data["kruskal_stat"] = np.nan
                    row_data["kruskal_pval"] = np.nan
                
                # ANOVA (parametric)
                valid_for_anova = [gv for gv in group_values.values() if len(gv) >= 2]
                if len(valid_for_anova) >= 2:
                    try:
                        stat, pval = stats.f_oneway(*valid_for_anova)
                        row_data["anova_stat"] = stat
                        row_data["anova_pval"] = pval
                    except Exception:
                        row_data["anova_stat"] = np.nan
                        row_data["anova_pval"] = np.nan
                else:
                    row_data["anova_stat"] = np.nan
                    row_data["anova_pval"] = np.nan
                
                # Welch ANOVA (parametric with unequal variances)
                if len(valid_for_anova) >= 2:
                    try:
                        # Welch ANOVA using scipy (one-way test with equal_var=False not available directly)
                        # We'll use a workaround: for 2 groups use Welch t-test, for 3+ use Kruskal as fallback
                        if len(groups) == 2:
                            g1_vals = group_values[groups[0]]
                            g2_vals = group_values[groups[1]]
                            if len(g1_vals) >= 2 and len(g2_vals) >= 2:
                                stat, pval = stats.ttest_ind(g1_vals, g2_vals, equal_var=False)
                                row_data["welch_stat"] = stat
                                row_data["welch_pval"] = pval
                            else:
                                row_data["welch_stat"] = np.nan
                                row_data["welch_pval"] = np.nan
                        else:
                            # For 3+ groups, Welch ANOVA is complex; use oneway with unequal var assumption
                            # scipy doesn't have direct Welch ANOVA, so we use Kruskal as robust alternative
                            row_data["welch_stat"] = row_data.get("kruskal_stat", np.nan)
                            row_data["welch_pval"] = row_data.get("kruskal_pval", np.nan)
                    except Exception:
                        row_data["welch_stat"] = np.nan
                        row_data["welch_pval"] = np.nan
                else:
                    row_data["welch_stat"] = np.nan
                    row_data["welch_pval"] = np.nan
                
                # Determine primary test based on choice
                if test_choice == "parametric":
                    row_data["primary_stat"] = row_data.get("anova_stat", np.nan)
                    row_data["primary_pval"] = row_data.get("anova_pval", np.nan)
                elif test_choice == "welch":
                    row_data["primary_stat"] = row_data.get("welch_stat", np.nan)
                    row_data["primary_pval"] = row_data.get("welch_pval", np.nan)
                else:  # nonparametric (default)
                    row_data["primary_stat"] = row_data.get("kruskal_stat", np.nan)
                    row_data["primary_pval"] = row_data.get("kruskal_pval", np.nan)
            else:
                row_data["kruskal_stat"] = np.nan
                row_data["kruskal_pval"] = np.nan
                row_data["anova_stat"] = np.nan
                row_data["anova_pval"] = np.nan
                row_data["welch_stat"] = np.nan
                row_data["welch_pval"] = np.nan
                row_data["primary_stat"] = np.nan
                row_data["primary_pval"] = np.nan
    
            # Pairwise tests
            for g1, g2 in combinations(groups, 2):
                v1, v2 = group_values.get(g1, []), group_values.get(g2, [])
                pair_key = f"{g1}_vs_{g2}"
                
                if len(v1) >= 1 and len(v2) >= 1:
                    # Mann-Whitney U (nonparametric, always compute)
                    try:
                        stat, pval = stats.mannwhitneyu(v1, v2, alternative="two-sided")
                        row_data[f"mwu_stat_{pair_key}"] = stat
                        row_data[f"mwu_pval_{pair_key}"] = pval
                    except Exception:
                        row_data[f"mwu_stat_{pair_key}"] = np.nan
                        row_data[f"mwu_pval_{pair_key}"] = np.nan
                    
                    # Student t-test (parametric, equal variances)
                    if len(v1) >= 2 and len(v2) >= 2:
                        try:
                            stat, pval = stats.ttest_ind(v1, v2, equal_var=True)
                            row_data[f"student_stat_{pair_key}"] = stat
                            row_data[f"student_pval_{pair_key}"] = pval
                        except Exception:
                            row_data[f"student_stat_{pair_key}"] = np.nan
                            row_data[f"student_pval_{pair_key}"] = np.nan
                    else:
                        row_data[f"student_stat_{pair_key}"] = np.nan
                        row_data[f"student_pval_{pair_key}"] = np.nan
                    
                    # Welch t-test (parametric, unequal variances)
                    if len(v1) >= 2 and len(v2) >= 2:
                        try:
                            stat, pval = stats.ttest_ind(v1, v2, equal_var=False)
                            row_data[f"welch_stat_{pair_key}"] = stat
                            row_data[f"welch_pval_{pair_key}"] = pval
                        except Exception:
                            row_data[f"welch_stat_{pair_key}"] = np.nan
                            row_data[f"welch_pval_{pair_key}"] = np.nan
                    else:
                        row_data[f"welch_stat_{pair_key}"] = np.nan
                        row_data[f"welch_pval_{pair_key}"] = np.nan
                    
                    # Effect size (rank-biserial for Mann-Whitney)
                    n1, n2 = len(v1), len(v2)
                    if n1 * n2 > 0 and not np.isnan(row_data.get(f"mwu_stat_{pair_key}", np.nan)):
                        row_data[f"effect_size_{pair_key}"] = 1 - 2 * row_data[f"mwu_stat_{pair_key}"] / (n1 * n2)
                    
                    # Log2 fold change of means
                    m1, m2 = np.mean(v1), np.mean(v2)
                    if m2 != 0 and m1 != 0:
                        row_data[f"log2fc_{pair_key}"] = np.log2(m1 / m2) if m1 > 0 and m2 > 0 else np.nan
                    row_data[f"diff_{pair_key}"] = m1 - m2
                else:
                    row_data[f"mwu_stat_{pair_key}"] = np.nan
                    row_data[f"mwu_pval_{pair_key}"] = np.nan
                    row_data[f"student_stat_{pair_key}"] = np.nan
                    row_data[f"student_pval_{pair_key}"] = np.nan
                    row_data[f"welch_stat_{pair_key}"] = np.nan
                    row_data[f"welch_pval_{pair_key}"] = np.nan
    
            rows.append(row_data)

    # Convert rows to DataFrame
    res_df = pd.DataFrame(rows)
    
    # Multiple testing correction (FDR Benjamini-Hochberg) for ALL p-value columns
    for col in res_df.columns:
        if col.endswith("_pval"):
            pvals = res_df[col].values
            mask = ~np.isnan(pvals)
            if mask.sum() > 0:
                _, corrected, _, _ = multipletests(pvals[mask], method="fdr_bh")
                adj_col = col.replace("_pval", "_padj")
                res_df[adj_col] = np.nan
                res_df.loc[mask, adj_col] = corrected

    res_df.to_csv(os.path.join(outdir, "differential_results.csv"), index=False)
    results["table"] = os.path.join(outdir, "differential_results.csv")
    results["stat_test_method"] = stat_test

    # Volcano-like plot for each pairwise comparison (using primary test)
    for g1, g2 in combinations(groups, 2):
        pair_key = f"{g1}_vs_{g2}"
        diff_col = f"diff_{pair_key}"
        
        # Choose p-value column based on test method
        if stat_test == "parametric":
            pval_col = f"student_pval_{pair_key}"
            padj_col = f"student_padj_{pair_key}"
            test_label = "Student"
        elif stat_test == "welch":
            pval_col = f"welch_pval_{pair_key}"
            padj_col = f"welch_padj_{pair_key}"
            test_label = "Welch"
        else:  # nonparametric or auto (default to nonparametric)
            pval_col = f"mwu_pval_{pair_key}"
            padj_col = f"mwu_padj_{pair_key}"
            test_label = "Mann-Whitney U"
        if diff_col in res_df.columns and padj_col in res_df.columns:
            pdf = res_df[[diff_col, padj_col]].dropna()
            if len(pdf) > 0:
                fig, ax = plt.subplots(figsize=(7, 5))
                neg_log_p = -np.log10(pdf[padj_col].clip(lower=1e-300))
                colors = ["red" if p < 0.05 else "grey" for p in pdf[padj_col]]
                ax.scatter(pdf[diff_col], neg_log_p, c=colors, alpha=0.6, s=20)
                ax.axhline(-np.log10(0.05), color="blue", linestyle="--", alpha=0.5)
                ax.set_xlabel(f"Difference ({g1} - {g2})")
                ax.set_ylabel("-log10(adjusted p-value)")
                ax.set_title(f"Volcano plot: {g1} vs {g2} ({test_label})")
                save_fig(fig, os.path.join(outdir, f"volcano_{pair_key}.png"))

    return results


# ---------------------------------------------------------------------------
# Multivariate Analysis (PCA, hierarchical clustering)
# ---------------------------------------------------------------------------

def multivariate_analysis(df, sample_cols, sample_info, outdir):
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # Transpose: samples as rows, BMKs as columns
    mat = ndf[sample_cols].T.copy()
    mat.columns = df["ID"].values if "ID" in df.columns else range(len(df))
    mat = mat.dropna(axis=1, how="all").fillna(0)

    if mat.shape[1] < 2 or mat.shape[0] < 2:
        log.warning("Not enough data for multivariate analysis")
        return results

    # PCA
    scaler = StandardScaler()
    scaled = scaler.fit_transform(mat)
    n_components = min(mat.shape[0], mat.shape[1], 10)
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(scaled)
    var_exp = pca.explained_variance_ratio_

    labels = [s["group"] for s in sample_info if s["col"] in mat.index]
    sample_labels = [s["sample"] for s in sample_info if s["col"] in mat.index]

    pca_df = pd.DataFrame(coords[:, :min(5, n_components)],
                          columns=[f"PC{i+1}" for i in range(min(5, n_components))],
                          index=mat.index)
    pca_df["group"] = labels
    pca_df["sample"] = sample_labels
    pca_df.to_csv(os.path.join(outdir, "pca_coordinates.csv"))

    # Variance explained
    var_df = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(var_exp))],
                           "variance_explained": var_exp,
                           "cumulative": np.cumsum(var_exp)})
    var_df.to_csv(os.path.join(outdir, "pca_variance.csv"), index=False)

    # Scree plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(1, len(var_exp) + 1), var_exp * 100, alpha=0.7, label="Individual")
    ax.plot(range(1, len(var_exp) + 1), np.cumsum(var_exp) * 100, "ro-", label="Cumulative")
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("PCA Scree Plot")
    ax.legend()
    save_fig(fig, os.path.join(outdir, "pca_scree.png"))

    # PCA biplot PC1 vs PC2
    if n_components >= 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        unique_groups = sorted(set(labels))
        colors = plt.cm.Set1(np.linspace(0, 1, max(len(unique_groups), 1)))
        for i, g in enumerate(unique_groups):
            mask = [l == g for l in labels]
            ax.scatter(coords[mask, 0], coords[mask, 1], c=[colors[i]], label=g, s=80, alpha=0.8)
            for j, m in enumerate(mask):
                if m:
                    ax.annotate(sample_labels[j], (coords[j, 0], coords[j, 1]),
                                fontsize=7, alpha=0.7)
        ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
        ax.set_title("PCA - Samples")
        ax.legend()
        save_fig(fig, os.path.join(outdir, "pca_biplot.png"))

    # Hierarchical clustering on samples
    if mat.shape[0] >= 2:
        try:
            dist = pdist(scaled, metric="euclidean")
            Z = linkage(dist, method="ward")
            fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.8), 5))
            dendrogram(Z, labels=[s.split("::")[1] for s in mat.index], ax=ax, leaf_rotation=90)
            ax.set_title("Hierarchical Clustering of Samples")
            save_fig(fig, os.path.join(outdir, "dendrogram_samples.png"))
        except Exception as e:
            log.warning(f"Dendrogram failed: {e}")

    # Sample correlation heatmap
    corr = mat.T.corr()
    corr.to_csv(os.path.join(outdir, "sample_correlation.csv"))
    fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.8), max(5, len(sample_cols) * 0.6)))
    short_labels = [s.split("::")[0] + "::" + s.split("::")[1] for s in corr.index]
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0, ax=ax,
                xticklabels=short_labels, yticklabels=short_labels)
    ax.set_title("Sample Correlation Matrix")
    save_fig(fig, os.path.join(outdir, "correlation_heatmap.png"))

    results["pca_variance"] = var_df.to_dict(orient="records")
    return results


# ---------------------------------------------------------------------------
# Correlation / Network Analysis
# ---------------------------------------------------------------------------

def correlation_network(df, sample_cols, outdir, max_bmks=50):
    """Compute BMK-BMK correlation matrix and heatmap.
    
    MEMORY FIX: Reduced max_bmks from 200 to 50 to prevent OOM.
    - Correlation matrix: N×N float64 = N² × 8 bytes
    - With N=200: 320 KB per section × 1491 sections = 477 MB accumulated
    - With N=50: 20 KB per section × 1491 sections = 30 MB accumulated (16x reduction)
    - Heatmap figure with N=50: ~4 MB instead of 64 MB per image
    """
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # BMK-BMK correlation (limit to top variable BMKs)
    mat = ndf[sample_cols].copy()
    mat.index = df["ID"].values if "ID" in df.columns else range(len(df))
    mat = mat.dropna(axis=0, how="all")
    
    # Early exit if too few BMKs
    if len(mat) < 3:
        log.info(f"  Skipping correlation analysis: only {len(mat)} BMKs (need at least 3)")
        return results
    
    var = mat.var(axis=1)
    top_idx = var.nlargest(min(max_bmks, len(var))).index
    sub = mat.loc[top_idx]

    if len(sub) >= 3:
        # Compute correlation matrix (memory intensive: N×N)
        corr = sub.T.corr()
        corr.to_csv(os.path.join(outdir, "bmk_correlation.csv"))

        # Limit figure size to prevent excessive memory usage
        fig_width = min(15, max(6, len(sub) * 0.15))   # Cap at 15 inches
        fig_height = min(12, max(5, len(sub) * 0.12))  # Cap at 12 inches
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(corr, cmap="RdBu_r", center=0, ax=ax,
                    xticklabels=len(sub) < 50, yticklabels=len(sub) < 50)
        ax.set_title(f"BMK Correlation (top {len(sub)} by variance)")
        save_fig(fig, os.path.join(outdir, "bmk_correlation_heatmap.png"))
        results["n_bmks_corr"] = len(sub)

    return results


# ---------------------------------------------------------------------------
# Feature Selection / Biomarker Ranking
# ---------------------------------------------------------------------------

def feature_ranking(df, sample_cols, sample_info, outdir, max_bmks=500, bmk_filter_cols=None):
    """Feature ranking using RandomForest importance.
    
    Args:
        bmk_filter_cols: list of column names to try in order for filtering significant BMKs
    """
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    groups = sorted(set(s["group"] for s in sample_info))
    results = {}

    if len(groups) < 2:
        return results

    # Build matrix: rows = samples, cols = BMKs
    mat = ndf[sample_cols].T.copy()
    bmk_ids = df["ID"].values if "ID" in df.columns else [str(i) for i in range(len(df))]
    mat.columns = bmk_ids
    mat = mat.fillna(0)

    labels = [group_for_col(sample_info, c) for c in sample_cols]
    le = LabelEncoder()
    y = le.fit_transform(labels)
    
    # Default filter column priority if not specified
    if bmk_filter_cols is None:
        bmk_filter_cols = ["primary_padj", "kruskal_padj", "welch_padj", "anova_padj"]
    
    # Filter BMKs: cascade through filter columns to collect significant ones
    # PASS 1: adjusted p-values (*_padj < 0.05)
    # PASS 2: raw p-values (*_pval < 0.05) if not enough
    diff_file = os.path.join(os.path.dirname(outdir), "differential", "differential_results.csv")
    selected_bmks = set()
    n_bmks_total = len(mat.columns)
    filter_log = []
    
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        
        # PASS 1: Cascade through adjusted p-values (strong evidence)
        for filter_col in bmk_filter_cols:
            if len(selected_bmks) >= max_bmks:
                break
                
            if filter_col in diff_df.columns and "ID" in diff_df.columns:
                sig_bmks = diff_df[diff_df[filter_col] < 0.05]["ID"].values
                valid_sig_bmks = [b for b in sig_bmks if b in mat.columns and b not in selected_bmks]
                remaining = max_bmks - len(selected_bmks)
                to_add = valid_sig_bmks[:remaining]
                
                if to_add:
                    selected_bmks.update(to_add)
                    filter_log.append(f"{filter_col}: +{len(to_add)}")
        
        # PASS 2: If not enough, cascade through raw p-values (suggestive evidence)
        if len(selected_bmks) < max_bmks:
            filter_log.append("|")
            for filter_col in bmk_filter_cols:
                if len(selected_bmks) >= max_bmks:
                    break
                
                # Convert *_padj to *_pval
                pval_col = filter_col.replace("_padj", "_pval")
                if pval_col != filter_col and pval_col in diff_df.columns and "ID" in diff_df.columns:
                    sig_bmks = diff_df[diff_df[pval_col] < 0.05]["ID"].values
                    valid_sig_bmks = [b for b in sig_bmks if b in mat.columns and b not in selected_bmks]
                    remaining = max_bmks - len(selected_bmks)
                    to_add = valid_sig_bmks[:remaining]
                    
                    if to_add:
                        selected_bmks.update(to_add)
                        filter_log.append(f"{pval_col}: +{len(to_add)}")
        
        if selected_bmks:
            selected_bmks = list(selected_bmks)
            log.info(f"  Collected {len(selected_bmks)} significant BMKs from {n_bmks_total} total for ranking")
            log.info(f"    Cascade: {' → '.join(filter_log)}")
        else:
            selected_bmks = None
    
    # If not enough significant BMKs, complete with top variable ones
    if selected_bmks is None or len(selected_bmks) < 10:
        if selected_bmks is None:
            selected_bmks = []
        var_scores = mat.var(axis=0)
        # Exclude already selected
        var_scores = var_scores[[b for b in var_scores.index if b not in selected_bmks]]
        remaining = max_bmks - len(selected_bmks)
        top_var = var_scores.nlargest(min(remaining, len(var_scores))).index.tolist()
        selected_bmks.extend(top_var)
        log.info(f"  Completed with {len(top_var)} top variable BMKs (total: {len(selected_bmks)} from {n_bmks_total})")
    
    # Filter matrix to selected BMKs
    mat_filtered = mat[selected_bmks].copy()

    # Random Forest importance
    if len(set(y)) >= 2 and mat_filtered.shape[0] >= 4 and mat_filtered.shape[1] >= 2:
        try:
            rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1, max_depth=10)
            rf.fit(mat_filtered, y)
            importances = rf.feature_importances_
            imp_df = pd.DataFrame({"bmk": mat_filtered.columns, "importance": importances})
            imp_df = imp_df.sort_values("importance", ascending=False)
            imp_df.to_csv(os.path.join(outdir, "rf_importance.csv"), index=False)
            results["rf_top10"] = imp_df.head(10).to_dict(orient="records")

            # Plot top 30
            top_n = min(30, len(imp_df))
            top = imp_df.head(top_n)
            fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.3)))
            ax.barh(range(top_n), top["importance"].values[::-1])
            ax.set_yticks(range(top_n))
            ax.set_yticklabels(top["bmk"].values[::-1], fontsize=7)
            ax.set_xlabel("Importance")
            ax.set_title(f"Random Forest Feature Importance (top 30 from {len(selected_bmks)} BMKs)")
            save_fig(fig, os.path.join(outdir, "rf_importance.png"))
        except Exception as e:
            log.warning(f"RF importance failed: {e}")

    # Variance-based ranking
    var_scores = ndf[sample_cols].var(axis=1)
    var_df = pd.DataFrame({"bmk": bmk_ids, "variance": var_scores.values})
    var_df = var_df.sort_values("variance", ascending=False)
    var_df.to_csv(os.path.join(outdir, "variance_ranking.csv"), index=False)

    # Kruskal-Wallis based ranking (from differential analysis if available)
    diff_file = os.path.join(os.path.dirname(outdir), "differential", "differential_results.csv")
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        if "kruskal_padj" in diff_df.columns:
            rank_df = diff_df[["ID", "kruskal_padj"]].dropna().sort_values("kruskal_padj")
            rank_df.to_csv(os.path.join(outdir, "kruskal_ranking.csv"), index=False)
            results["kruskal_top10"] = rank_df.head(10).to_dict(orient="records")

    return results


# ---------------------------------------------------------------------------
# Classification / Predictive Modeling
# ---------------------------------------------------------------------------

def classification_analysis(df, sample_cols, sample_info, outdir, max_bmks=500, bmk_filter_cols=None):
    """Classification analysis with multiple classifiers.
    
    Args:
        bmk_filter_cols: list of column names to try in order for filtering significant BMKs
    """
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    groups = sorted(set(s["group"] for s in sample_info))
    results = {}

    if len(groups) < 2:
        return results

    mat = ndf[sample_cols].T.fillna(0).copy()
    bmk_ids = df["ID"].values if "ID" in df.columns else [str(i) for i in range(len(df))]
    mat.columns = bmk_ids
    labels = [group_for_col(sample_info, c) for c in sample_cols]
    le = LabelEncoder()
    y = le.fit_transform(labels)

    n_samples = len(y)
    n_classes = len(set(y))
    n_bmks_total = mat.shape[1]

    if n_samples < 4 or n_classes < 2:
        log.warning("Not enough samples for classification")
        return results
    
    # Default filter column priority if not specified
    if bmk_filter_cols is None:
        bmk_filter_cols = ["primary_padj", "kruskal_padj", "welch_padj", "anova_padj"]
    
    # Filter BMKs: cascade through filter columns to collect significant ones
    diff_file = os.path.join(os.path.dirname(outdir), "differential", "differential_results.csv")
    selected_bmks = set()
    filter_log = []
    
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        
        # Cascade through filter columns in priority order
        for filter_col in bmk_filter_cols:
            if len(selected_bmks) >= max_bmks:
                break  # Stop if we have enough
                
            if filter_col in diff_df.columns and "ID" in diff_df.columns:
                # Get significant BMKs from this column
                sig_bmks = diff_df[diff_df[filter_col] < 0.05]["ID"].values
                valid_sig_bmks = [b for b in sig_bmks if b in mat.columns and b not in selected_bmks]
                
                # Add up to remaining slots
                remaining = max_bmks - len(selected_bmks)
                to_add = valid_sig_bmks[:remaining]
                
                if to_add:
                    selected_bmks.update(to_add)
                    filter_log.append(f"{filter_col}: +{len(to_add)} BMKs")
        
        if selected_bmks:
            selected_bmks = list(selected_bmks)
            log.info(f"  Collected {len(selected_bmks)} significant BMKs from {n_bmks_total} total for classification")
            log.info(f"    Cascade: {' → '.join(filter_log)}")
        else:
            selected_bmks = None
    
    # If not enough significant BMKs, complete with top variable ones
    if selected_bmks is None or len(selected_bmks) < 10:
        if selected_bmks is None:
            selected_bmks = []
        var_scores = mat.var(axis=0)
        # Exclude already selected
        var_scores = var_scores[[b for b in var_scores.index if b not in selected_bmks]]
        remaining = max_bmks - len(selected_bmks)
        top_var = var_scores.nlargest(min(remaining, len(var_scores))).index.tolist()
        selected_bmks.extend(top_var)
        log.info(f"  Completed with {len(top_var)} top variable BMKs (total: {len(selected_bmks)} from {n_bmks_total})")
    
    # Filter matrix to selected BMKs
    mat = mat[selected_bmks].copy()
    n_bmks = mat.shape[1]
    
    if n_bmks < max(2, n_classes):
        log.warning(f"Not enough biomarkers for classification ({n_bmks} BMKs after filtering, need at least {max(2, n_classes)})")
        return results

    # Use LOO or stratified k-fold depending on sample count
    if n_samples < 10:
        cv = LeaveOneOut()
        cv_name = "LOO"
    else:
        cv = StratifiedKFold(n_splits=min(5, n_samples), shuffle=True, random_state=42)
        cv_name = "StratifiedKFold"

    classifiers = {}
    classifiers["RandomForest"] = RandomForestClassifier(n_estimators=50, random_state=42, n_jobs=1, max_depth=10)

    if n_classes == 2:
        classifiers["GradientBoosting"] = GradientBoostingClassifier(n_estimators=50, random_state=42)

    # LDA only if feasible
    if n_samples > n_classes and mat.shape[1] > 0:
        try:
            classifiers["LDA"] = LinearDiscriminantAnalysis()
        except Exception:
            pass

    clf_results = {}
    for name, clf in classifiers.items():
        try:
            scores = cross_val_score(clf, mat, y, cv=cv, scoring="accuracy")
            clf_results[name] = {
                "mean_accuracy": float(np.mean(scores)),
                "std_accuracy": float(np.std(scores)),
                "cv_method": cv_name,
                "scores": scores.tolist(),
                "n_bmks_used": n_bmks,
            }
        except Exception as e:
            log.warning(f"Classifier {name} failed: {e}")

    results["classifiers"] = clf_results
    clf_df = pd.DataFrame([
        {"classifier": k, "mean_accuracy": v["mean_accuracy"],
         "std_accuracy": v["std_accuracy"], "cv_method": v["cv_method"],
         "n_bmks_used": v.get("n_bmks_used", n_bmks)}
        for k, v in clf_results.items()
    ])
    clf_df.to_csv(os.path.join(outdir, "classification_results.csv"), index=False)

    # Plot
    if clf_results:
        fig, ax = plt.subplots(figsize=(6, 4))
        names = list(clf_results.keys())
        means = [clf_results[n]["mean_accuracy"] for n in names]
        stds = [clf_results[n]["std_accuracy"] for n in names]
        ax.bar(names, means, yerr=stds, alpha=0.7, capsize=5)
        ax.set_ylabel("Accuracy")
        ax.set_title(f"Classification Accuracy ({cv_name}, {n_bmks} BMKs)")
        ax.set_ylim(0, 1.1)
        save_fig(fig, os.path.join(outdir, "classification_accuracy.png"))

    return results


# ---------------------------------------------------------------------------
# Stability / Robustness (replicate concordance)
# ---------------------------------------------------------------------------

def stability_analysis(df, sample_cols, sample_info, outdir):
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # Group by (group, sample) and check replicate correlation
    rep_groups = defaultdict(list)
    for s in sample_info:
        key = (s["group"], s["sample"])
        rep_groups[key].append(s["col"])

    rep_corrs = []
    for key, cols in rep_groups.items():
        if len(cols) >= 2:
            for c1, c2 in combinations(cols, 2):
                v1 = pd.to_numeric(ndf[c1], errors="coerce")
                v2 = pd.to_numeric(ndf[c2], errors="coerce")
                mask = v1.notna() & v2.notna()
                if mask.sum() >= 3:
                    r, p = stats.pearsonr(v1[mask], v2[mask])
                    rep_corrs.append({
                        "group": key[0], "sample": key[1],
                        "rep1": c1, "rep2": c2,
                        "pearson_r": r, "pval": p
                    })

    if rep_corrs:
        rc_df = pd.DataFrame(rep_corrs)
        rc_df.to_csv(os.path.join(outdir, "replicate_correlation.csv"), index=False)
        results["replicate_correlations"] = rc_df.to_dict(orient="records")

    # Coefficient of variation per BMK across all samples
    cv_vals = ndf[sample_cols].apply(lambda row: row.std() / row.mean() if row.mean() != 0 else np.nan, axis=1)
    cv_df = pd.DataFrame({"bmk": df["ID"].values if "ID" in df.columns else range(len(df)),
                           "cv": cv_vals.values})
    cv_df = cv_df.sort_values("cv")
    cv_df.to_csv(os.path.join(outdir, "coefficient_variation.csv"), index=False)

    # CV histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(cv_df["cv"].dropna(), bins=30, alpha=0.7, edgecolor="black")
    ax.set_xlabel("Coefficient of Variation")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of CV across BMKs")
    save_fig(fig, os.path.join(outdir, "cv_histogram.png"))

    return results


# ---------------------------------------------------------------------------
# Batch Effect Detection
# ---------------------------------------------------------------------------

def batch_effect_analysis(df, sample_cols, sample_info, outdir):
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    results = {}

    # Check if samples cluster by batch (sample) rather than group
    # Use PCA and check if samples separate by sample-id vs group
    mat = ndf[sample_cols].T.fillna(0)
    if mat.shape[0] < 3 or mat.shape[1] < 2:
        return results

    scaler = StandardScaler()
    scaled = scaler.fit_transform(mat)
    pca = PCA(n_components=min(3, mat.shape[0], mat.shape[1]))
    coords = pca.fit_transform(scaled)

    groups = [group_for_col(sample_info, c) for c in sample_cols]
    samples = [s["sample"] for s in sample_info if s["col"] in sample_cols[:len(groups)]]

    # Plot colored by sample (batch)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    unique_groups = sorted(set(groups))
    colors_g = plt.cm.Set1(np.linspace(0, 1, max(len(unique_groups), 1)))
    for i, g in enumerate(unique_groups):
        mask = [l == g for l in groups]
        axes[0].scatter(coords[mask, 0], coords[mask, 1], c=[colors_g[i]], label=g, s=80)
    axes[0].set_title("PCA colored by Group")
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].legend()

    unique_samples = sorted(set(samples))
    colors_s = plt.cm.Set2(np.linspace(0, 1, max(len(unique_samples), 1)))
    for i, s in enumerate(unique_samples):
        mask = [l == s for l in samples]
        axes[1].scatter(coords[mask, 0], coords[mask, 1], c=[colors_s[i]], label=s, s=80)
    axes[1].set_title("PCA colored by Sample (batch)")
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")
    axes[1].legend(fontsize=7)

    save_fig(fig, os.path.join(outdir, "batch_effect_pca.png"))
    return results


# ---------------------------------------------------------------------------
# Comprehensive heatmap for a section
# ---------------------------------------------------------------------------

def section_heatmap(df, sample_cols, sample_info, outdir, title="", max_rows=100):
    safe_mkdir(outdir)
    ndf = numeric_df(df, sample_cols)
    
    # Sort sample_cols by group, then rep
    sample_order = sorted(sample_info, key=lambda s: (s["group"], s["rep"], s["sample"]))
    sorted_cols = [s["col"] for s in sample_order if s["col"] in sample_cols]
    
    mat = ndf[sorted_cols].copy()
    mat = mat.dropna(how="all")
    
    # CRITICAL: Limit rows BEFORE any processing to avoid matplotlib memory errors
    # With 200k+ BMKs, matplotlib cannot create images (pixel limit: 2^23 per dimension)
    if len(mat) > max_rows:
        var = mat.var(axis=1)
        mat = mat.loc[var.nlargest(max_rows).index]
    
    if len(mat) == 0:
        return
    
    # Create descriptive Y-axis labels
    # Add SeqID/Ptype/Mode affixes when multiple values exist
    if "Mode" in df.columns and "Ctype" in df.columns:
        # Check if we need SeqID/Ptype/Mode affixes (mixed BMK types or chromosomes)
        df_subset = df.loc[mat.index]
        has_multiple_seqids = "SeqID" in df_subset.columns and df_subset["SeqID"].nunique() > 1
        has_multiple_ptypes = "Ptype" in df_subset.columns and df_subset["Ptype"].nunique() > 1
        has_multiple_modes = "Mode" in df_subset.columns and df_subset["Mode"].nunique() > 1
        add_affixes = has_multiple_seqids or has_multiple_ptypes or has_multiple_modes
        
        ylabels = []
        sort_keys = []  # For custom sorting: (priority, label)
        for idx in mat.index:
            mode = df.loc[idx, "Mode"] if "Mode" in df.columns else ""
            ctype = df.loc[idx, "Ctype"] if "Ctype" in df.columns else ""
            ptype = df.loc[idx, "Ptype"] if "Ptype" in df.columns else ""
            seqid = df.loc[idx, "SeqID"] if "SeqID" in df.columns else ""
            
            # Build base label
            if mode == "all_sites":
                base_label = "all_sites"
            elif ctype and ctype != ".":
                base_label = str(ctype)
            elif "ID" in df.columns:
                base_label = str(df.loc[idx, "ID"])
            else:
                base_label = str(idx)
            
            # Add affixes only if values vary across rows
            if add_affixes:
                # SeqID as PREFIX (first, if multiple chromosomes)
                if has_multiple_seqids and seqid and seqid != ".":
                    label = f"{seqid}_{base_label}"
                else:
                    label = base_label
                
                # Ptype as PREFIX or SUFFIX depending on context
                if ptype and ptype != ".":
                    # If Ptypes vary OR SeqIDs vary, add Ptype
                    if has_multiple_ptypes or has_multiple_seqids:
                        if has_multiple_seqids:
                            # For all_sequence_bmks: SeqID_Ptype_Ctype
                            label = f"{seqid}_{ptype}_{base_label}" if seqid and seqid != "." else f"{ptype}_{base_label}"
                        else:
                            # For all_global_bmks: Ptype_Ctype
                            label = f"{ptype}_{base_label}"
                
                # Mode as SUFFIX (after), if not "all_sites"
                if mode and mode != "all_sites":
                    label += f"_{mode}"
                
                # Determine sort priority for mixed BMKs:
                # 1. all_sites first
                # 2. no Ptype (Ptype == ".")
                # 3. rest alphabetically
                if mode == "all_sites":
                    priority = 0
                elif not ptype or ptype == ".":
                    priority = 1
                else:
                    priority = 2
                sort_keys.append((priority, label, idx))
            else:
                label = base_label
                sort_keys.append((0, label, idx))
            
            ylabels.append(label)
        
        # Sort rows if affixes were added (mixed BMK types)
        if add_affixes and len(sort_keys) > 1:
            # Sort by (priority, label alphabetically)
            sorted_items = sorted(sort_keys, key=lambda x: (x[0], x[1]))
            sorted_indices = [item[2] for item in sorted_items]
            mat = mat.loc[sorted_indices]
            ylabels = [item[1] for item in sorted_items]
        
        mat.index = ylabels

    mat = mat.dropna(how="all")
    if len(mat) == 0:
        return

    if len(mat) > max_rows:
        # Keep most variable
        var = mat.var(axis=1)
        mat = mat.loc[var.nlargest(max_rows).index]

    fig_h = max(4, len(mat) * 0.25)
    fig_w = max(6, len(sorted_cols) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    # Create labels: group::sample::rep for x-axis
    xlabels = [f"{s['group']}::{s['sample']}::{s['rep']}" for s in sample_order if s['col'] in sorted_cols]
    
    # Determine y-axis label display and font size
    show_ylabels = len(mat) < 200  # Show labels up to 200 rows
    yticklabel_fontsize = max(4, min(8, 500 / len(mat)))  # Adaptive font size
    
    sns.heatmap(mat.astype(float), cmap="YlOrRd", ax=ax,
                xticklabels=xlabels, 
                yticklabels=show_ylabels,
                linewidths=0.5 if len(mat) < 30 else 0,
                cbar_kws={'label': 'Value'})
    
    if show_ylabels:
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=yticklabel_fontsize, rotation=0)
    
    ax.set_xlabel("Sample (Group::Sample::Replicate)", fontsize=10)
    ax.set_ylabel("Biomarker Type" if show_ylabels else f"Biomarker ({len(mat)} total)", fontsize=10)
    ax.set_title(title or "Heatmap")
    plt.xticks(rotation=90, fontsize=8)
    save_fig(fig, os.path.join(outdir, "heatmap.png"))


# ---------------------------------------------------------------------------
# Top-level orchestration per section
# ---------------------------------------------------------------------------

def analyze_section(df, sample_cols, sample_info, outdir, section_name, stat_test="auto", bmk_filter_cols=None, max_bmks=500):
    """Run all analyses on a filtered dataframe for a given report section."""
    len_df = len(df)
    len_sample_cols = len(sample_cols)
    if len_df == 0:
        log.info(f"  Section '{section_name}' is empty, skipping.")
        return {}

    # Create log prefix for this section
    log_prefix = f"[{section_name[:50]}]"  # Limit to 50 chars
    
    log.info(f"{log_prefix} Analyzing section: {section_name} ({len_df} BMKs, {len_sample_cols} samples)")
    safe_mkdir(outdir)
    results = {"n_bmks": len_df, "n_samples": len_sample_cols, "section": section_name}
    
    # EARLY PRE-FILTER: For very large sections (>3k BMKs), reduce to top 3k by variance
    # This prevents memory crashes in parallel workers for all_sequence_bmks / all_global_bmks
    MAX_BMKS_FOR_ANALYSIS = 3000  # Reduced from 5000 to prevent hangs
    if len_df > MAX_BMKS_FOR_ANALYSIS:
        log.info(f"{log_prefix} {len_df} BMKs exceeds {MAX_BMKS_FOR_ANALYSIS} limit")
        log.info(f"{log_prefix} Pre-filtering to top {MAX_BMKS_FOR_ANALYSIS} BMKs by variance to prevent memory issues...")
        ndf_prefilter = numeric_df(df, sample_cols)
        variance_prefilter = ndf_prefilter[sample_cols].var(axis=1)
        top_indices = variance_prefilter.nlargest(MAX_BMKS_FOR_ANALYSIS).index
        df = df.loc[top_indices].reset_index(drop=True).copy()
        len_df = len(df)
        log.info(f"{log_prefix} Pre-filtered to {len_df} BMKs")
        results["n_bmks_prefiltered"] = len_df

    # Save filtered data
    df.to_csv(os.path.join(outdir, "data.csv"), index=False)

    # 1. QC
    log.info(f"{log_prefix} QC analysis...")
    try:
        results["qc"] = qc_analysis(df.reset_index(drop=True), sample_cols, os.path.join(outdir, "qc"))
    except Exception as e:
        log.warning(f"{log_prefix} QC failed: {e}")

    # 2. Descriptive stats
    log.info(f"{log_prefix} Descriptive stats...")
    try:
        results["descriptive"] = descriptive_stats(df.reset_index(drop=True), sample_cols, sample_info, os.path.join(outdir, "descriptive"))
    except Exception as e:
        log.warning(f"{log_prefix} Descriptive stats failed: {e}")

    # PRE-FILTER: Remove BMKs with near-zero variance for all subsequent analyses
    # This prevents numerical issues in: differential tests, PCA, correlation, heatmaps, ML models
    # QC and descriptive stats above use ALL BMKs to identify constants
    ndf_tmp = numeric_df(df, sample_cols)
    variance = ndf_tmp[sample_cols].var(axis=1)
    variable_mask = variance > 1e-10
    n_filtered_out = (~variable_mask).sum()
    
    if n_filtered_out > 0:
        log.info(f"{log_prefix} Removing {n_filtered_out} BMKs with near-zero variance (keeping {variable_mask.sum()} variable)")
        df_variable = df[variable_mask].reset_index(drop=True).copy()
        n_bmks_variable = len(df_variable)
    else:
        df_variable = df.copy()
        n_bmks_variable = len_df
    
    # Early exit if no variable BMKs
    if n_bmks_variable == 0:
        log.warning(f"{log_prefix} No variable BMKs remaining after pre-filtering")
        return results
    
    # Update results with variable BMK count
    results["n_bmks_variable"] = n_bmks_variable

    # 3. Differential analysis (use variable BMKs only)
    log.info(f"{log_prefix} Differential analysis ({n_bmks_variable} BMKs)...")
    try:
        results["differential"] = differential_analysis(df_variable, sample_cols, sample_info, os.path.join(outdir, "differential"), stat_test=stat_test)
    except Exception as e:
        log.warning(f"{log_prefix} Differential analysis failed: {e}")

    # 4. Multivariate (use variable BMKs only)
    log.info(f"{log_prefix} Multivariate analysis ({n_bmks_variable} BMKs)...")
    try:
        results["multivariate"] = multivariate_analysis(df_variable, sample_cols, sample_info, os.path.join(outdir, "multivariate"))
    except Exception as e:
        log.warning(f"{log_prefix} Multivariate analysis failed: {e}")

    # 5. Correlation / Network (use variable BMKs only)
    log.info(f"{log_prefix} Correlation / Network analysis ({n_bmks_variable} BMKs)...")
    try:
        results["correlation"] = correlation_network(df_variable, sample_cols, os.path.join(outdir, "correlation"))
    except Exception as e:
        log.warning(f"{log_prefix} Correlation analysis failed: {e}")

    # 6. Feature ranking (use variable BMKs only)
    log.info(f"{log_prefix} Feature ranking ({n_bmks_variable} BMKs)...")
    try:
        results["ranking"] = feature_ranking(df_variable, sample_cols, sample_info, os.path.join(outdir, "ranking"), max_bmks=max_bmks, bmk_filter_cols=bmk_filter_cols)
    except Exception as e:
        log.warning(f"{log_prefix} Feature ranking failed: {e}")

    # 7. Classification (use variable BMKs only)
    log.info(f"{log_prefix} Classification ({n_bmks_variable} BMKs)...")
    try:
        results["classification"] = classification_analysis(df_variable, sample_cols, sample_info, os.path.join(outdir, "classification"), max_bmks=max_bmks, bmk_filter_cols=bmk_filter_cols)
    except Exception as e:
        log.warning(f"{log_prefix} Classification failed: {e}")

    # 8. Stability (use variable BMKs only)
    log.info(f"{log_prefix} Stability analysis ({n_bmks_variable} BMKs)...")
    try:
        results["stability"] = stability_analysis(df_variable, sample_cols, sample_info, os.path.join(outdir, "stability"))
    except Exception as e:
        log.warning(f"{log_prefix} Stability analysis failed: {e}")

    # 9. Batch effect (use variable BMKs only)
    log.info(f"{log_prefix} Batch effect analysis ({n_bmks_variable} BMKs)...")
    try:
        results["batch"] = batch_effect_analysis(df_variable, sample_cols, sample_info, os.path.join(outdir, "batch"))
    except Exception as e:
        log.warning(f"{log_prefix} Batch effect analysis failed: {e}")

    # 10. Heatmap (use variable BMKs only)
    log.info(f"{log_prefix} Heatmap generation ({n_bmks_variable} BMKs)...")
    try:
        section_heatmap(df_variable, sample_cols, sample_info, outdir, title=section_name)
        log.info(f"{log_prefix} Heatmap generation completed")
    except Exception as e:
        log.warning(f"{log_prefix} Heatmap failed: {e}")

    log.info(f"{log_prefix} Section analysis completed successfully")
    return results


# ---------------------------------------------------------------------------
# Build feature hierarchy tree
# ---------------------------------------------------------------------------

def build_feature_tree(features_df):
    """Build parent-child tree for features.

    Each row's 'direct parent' is the last element in the comma-separated
    ParentIDs column (the first element is always '.').
    Top-level features have ParentIDs == '.'.
    """
    tree = {}  # id -> {"row": row, "children": []}
    for _, row in features_df.iterrows():
        fid = row["ID"]
        tree[fid] = {"data": row, "children": []}

    for _, row in features_df.iterrows():
        fid = row["ID"]
        parents = str(row["ParentIDs"])
        if parents == ".":
            continue  # top-level
        parts = [p.strip() for p in parents.split(",") if p.strip() != "."]
        if parts:
            direct_parent = parts[-1]
            if direct_parent in tree:
                tree[direct_parent]["children"].append(fid)

    # Find top-level features
    top_ids = [fid for fid, node in tree.items()
               if str(features_df.loc[features_df["ID"] == fid, "ParentIDs"].iloc[0]) == "."]

    return tree, top_ids


# ---------------------------------------------------------------------------
# Parallel execution wrapper
# ---------------------------------------------------------------------------

def analyze_section_wrapper(args_tuple):
    """Wrapper for analyze_section to be used with ProcessPoolExecutor.
    
    Args:
        args_tuple: (df_source, sample_cols, sample_info, outdir, section_name, result_key, stat_test, bmk_filter_cols, max_bmks)
            where df_source can be:
            - dict: legacy format, converted back to DataFrame
            - tuple: ('pickled', bytes_data) for memory-efficient transfer
        
    Returns:
        (result_key, results_dict)
    """
    import traceback
    import signal
    import gc
    import pickle
    import random
    import time as time_module
    df_source, sample_cols, sample_info, outdir, section_name, result_key, stat_test, bmk_filter_cols, max_bmks = args_tuple
    
    # OPTIMIZATION: Add startup jitter to avoid synchronized worker restarts.
    # With spawn + max_tasks_per_child=20, all workers start together and may
    # finish their 20th task at the same time → 6 workers restart simultaneously
    # → 6 × 275MB = 1.65GB spike. Jitter spreads restarts over 2 seconds.
    time_module.sleep(random.uniform(0, 2))
    
    # Log at start with worker PID and memory
    pid = os.getpid()
    mem_start = None
    if psutil:
        try:
            proc = psutil.Process(pid)
            mem_start = proc.memory_info().rss / (1024 * 1024)  # MB
            log.info(f"  ▶ START [PID {pid}, {mem_start:.0f}MB]: {section_name}")
        except:
            log.info(f"  ▶ START [PID {pid}]: {section_name}")
    else:
        log.info(f"  ▶ START [PID {pid}]: {section_name}")
    
    # Setup timeout alarm (120 seconds max per task)
    def timeout_handler(signum, frame):
        raise TimeoutError(f"Task exceeded 120 seconds: {section_name}")
    
    try:
        # Set alarm for 120 seconds
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(120)
        
        # Reconstruct DataFrame from source (dict or pickled bytes)
        if isinstance(df_source, tuple) and df_source[0] == 'pickled':
            # Memory-efficient: unpickle compressed data
            import pickle
            df = pickle.loads(df_source[1])
            log.info(f"  ▶ [PID {pid}] DataFrame unpickled: {len(df)} rows")
        else:
            # Legacy: dict format
            df = pd.DataFrame(df_source)
            log.info(f"  ▶ [PID {pid}] DataFrame reconstructed: {len(df)} rows")
        
        # Run analysis
        results = analyze_section(df, sample_cols, sample_info, outdir, section_name, stat_test=stat_test, bmk_filter_cols=bmk_filter_cols, max_bmks=max_bmks)
        
        # Cancel alarm
        signal.alarm(0)
        
        # MEMORY FIX #2: Return slim results dict (only file paths, not full data)
        # Full results dict is ~500KB per task × 1491 tasks = 745MB accumulated in all_results
        # global_ranking() only needs the differential.table path (50 bytes)
        slim_results = {
            "differential_table": os.path.join(outdir, "differential", "differential_results.csv")
        }
        
        # Explicitly clean up to help garbage collector (critical for parallel execution)
        del df
        del df_source
        del results  # MEMORY FIX #3: Immediate cleanup of full results dict in worker
        gc.collect()
        
        # Also force matplotlib to clean up any lingering figures
        plt.close('all')
        
        # Log completion with memory
        if psutil and mem_start:
            try:
                proc = psutil.Process(pid)
                mem_end = proc.memory_info().rss / (1024 * 1024)
                log.info(f"  ✓ DONE [PID {pid}, {mem_end:.0f}MB, Δ{mem_end-mem_start:+.0f}MB]: {section_name}")
            except:
                log.info(f"  ✓ DONE [PID {pid}]: {section_name}")
        else:
            log.info(f"  ✓ DONE [PID {pid}]: {section_name}")
        
        return result_key, slim_results
        
    except TimeoutError as e:
        signal.alarm(0)  # Cancel alarm
        log.error(f"  ✗ TIMEOUT [PID {pid}]: {section_name} - {e}")
        # Clean up even on error
        try:
            del df
            del df_source
            plt.close('all')
            gc.collect()
        except:
            pass
        raise
    except Exception as e:
        signal.alarm(0)  # Cancel alarm
        # Log full error details before crash
        log.error(f"  ✗ CRASH [PID {pid}]: {section_name}")
        log.error(f"     Error: {type(e).__name__}: {e}")
        log.error(f"     Stacktrace:\n{traceback.format_exc()}")
        # Clean up even on error
        try:
            del df
            del df_source
            plt.close('all')
            gc.collect()
        except:
            pass
        raise


# ---------------------------------------------------------------------------
# Global Biomarker Ranking (across all sections)
# ---------------------------------------------------------------------------

def global_ranking(all_results, outdir):
    """Aggregate rankings across sections to produce a global ranking."""
    safe_mkdir(outdir)

    # Collect differential results
    diff_files = []
    for vtype, vdata in all_results.items():
        for mtype, mdata in vdata.items():
            for section, sdata in mdata.items():
                # MEMORY FIX: sdata is now a slim dict with only "differential_table" key
                diff_path = sdata.get("differential_table")
                if diff_path and os.path.exists(diff_path):
                    ddf = pd.read_csv(diff_path)
                    # Defragment DataFrame before adding columns to avoid PerformanceWarning
                    ddf = ddf.copy()
                    # Add metadata columns
                    ddf["value_type"] = vtype
                    ddf["mtype"] = mtype
                    ddf["section"] = section
                    diff_files.append(ddf)

    if diff_files:
        all_diff = pd.concat(diff_files, ignore_index=True)
        # Rank by kruskal_padj
        if "kruskal_padj" in all_diff.columns:
            ranked = all_diff.dropna(subset=["kruskal_padj"]).sort_values("kruskal_padj")
            ranked.to_csv(os.path.join(outdir, "global_ranking_kruskal.csv"), index=False)

            # Top 50 plot
            top_n = min(50, len(ranked))
            top = ranked.head(top_n)
            fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.3)))
            neg_log_p = -np.log10(top["kruskal_padj"].clip(lower=1e-300))
            labels = top["ID"].astype(str) + " [" + top["section"].astype(str) + "]"
            ax.barh(range(top_n), neg_log_p.values[::-1])
            ax.set_yticks(range(top_n))
            ax.set_yticklabels(labels.values[::-1], fontsize=6)
            ax.set_xlabel("-log10(adjusted p-value)")
            ax.set_title("Global BMK Ranking by Significance (top 50)")
            
            # Add reference lines for p-value thresholds
            ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7, label='p=0.05')
            ax.axvline(-np.log10(0.1), color='orange', linestyle='--', linewidth=1, alpha=0.7, label='p=0.1')
            ax.legend(loc='lower right', fontsize=8)
            
            save_fig(fig, os.path.join(outdir, "global_ranking_plot.png"))

    return outdir


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="barometer Biomarker Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
STATISTICAL TESTS:
  The pipeline ALWAYS computes ALL differential tests for each biomarker:
  
  Global tests (comparing all groups):
    • Kruskal-Wallis → kruskal_pval, kruskal_padj
    • ANOVA         → anova_pval, anova_padj
    • Welch ANOVA   → welch_pval, welch_padj
  
  Pairwise tests (comparing 2 groups X vs Y):
    • Mann-Whitney U (= Wilcoxon rank-sum) → mwu_pval_X_vs_Y, mwu_padj_X_vs_Y
    • Student t-test                       → student_pval_X_vs_Y, student_padj_X_vs_Y
    • Welch t-test                         → welch_pval_X_vs_Y, welch_padj_X_vs_Y
  
  Assumption tests (only with --stat-test auto):
    • Shapiro-Wilk → shapiro_pvals (per group)
    • Bartlett     → bartlett_pval
  
  All tests receive FDR correction (Benjamini-Hochberg) → *_padj columns
  All results are saved in differential_results.csv
  
  OPTION 1: --stat-test (controls which test is emphasized in plots)
    This creates primary_pval and primary_padj columns that COPY the selected test:
    
    nonparametric : primary_* = kruskal_* (DEFAULT, robust)
    parametric    : primary_* = anova_*
    welch         : primary_* = welch_*
    auto          : primary_* = auto-selected test (varies per BMK)
    kruskal       : alias for nonparametric
  
  Example: --stat-test welch creates primary_padj as a copy of welch_padj
  Volcano plots use the primary_* columns for visualization.
  
  OPTION 2: --bmk-filter (cascade priority for selecting significant BMKs)
    Accepts a PRIORITY LIST of *_padj columns to collect significant BMKs
    (padj < 0.05) for RandomForest and classification.
    
    The algorithm works in TWO-PASS CASCADE:
      PASS 1 (adjusted p-values, strongest evidence):
        1. Take all BMKs with 1st column *_padj < 0.05
        2. If < --max-bmks, add BMKs from 2nd column *_padj < 0.05 (not already selected)
        3. Continue until reaching --max-bmks or exhausting all *_padj columns
      
      PASS 2 (raw p-values, suggestive evidence, only if PASS 1 insufficient):
        4. Convert columns to *_pval (e.g., kruskal_padj → kruskal_pval)
        5. Add BMKs with *_pval < 0.05 (not already selected)
        6. Continue until reaching --max-bmks
      
      PASS 3 (if still insufficient):
        7. Complete with top variable BMKs (by variance)
    
    Default cascade: primary_padj → kruskal_padj → welch_padj → anova_padj
    Then (if needed): primary_pval → kruskal_pval → welch_pval → anova_pval
    
    Available columns (choose your priority order):
      • primary_padj      : controlled by --stat-test (recommended)
      • kruskal_padj      : Kruskal-Wallis (robust, non-parametric)
      • anova_padj        : ANOVA (parametric, equal variances)
      • welch_padj        : Welch (parametric, unequal variances)
      • mwu_padj_X_vs_Y   : Mann-Whitney U for specific pair
      • student_padj_X_vs_Y : Student t-test for specific pair
    
    Examples:
      --bmk-filter kruskal_padj              : only use Kruskal (padj then pval)
      --bmk-filter primary_padj kruskal_padj : cascade from primary to kruskal
      --bmk-filter welch_padj anova_padj     : prioritize Welch, fallback to ANOVA
  
  OPTION 3: --max-bmks (maximum BMKs for ML models)
    Sets the hard limit for RandomForest and classification (default: 500).
    Prevents memory crashes with large datasets (e.g., 40,000+ BMKs).
    Higher values = more accurate but slower and more memory-intensive.

EXAMPLES:
  # Default: cascade adjusted + raw p-values, max 500 BMKs:
  ./barometer_analyze.py -a data.tsv -o results/ -j 4
  
  # Use only Kruskal-Wallis (padj then pval if needed):
  ./barometer_analyze.py -a data.tsv -o results/ --bmk-filter kruskal_padj
  
  # Prioritize Welch, fallback to Kruskal, use 1000 BMKs:
  ./barometer_analyze.py -a data.tsv -o results/ \
    --stat-test welch --bmk-filter welch_padj kruskal_padj --max-bmks 1000
  
  # Conservative: only most significant BMKs (200 max):
  ./barometer_analyze.py -a data.tsv -o results/ --max-bmks 200
  
  # Aggressive: maximize BMK collection with high limit:
  ./barometer_analyze.py -a data.tsv -o results/ --max-bmks 1500
        """)
    parser.add_argument("-a", "--aggregates", default=None, help="Aggregates TSV file (optional)")
    parser.add_argument("-f", "--features", default=None, help="Features TSV file (optional)")
    parser.add_argument("-o", "--outdir", default="barometer_results", help="Output directory")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Number of parallel jobs (default: 1, use -1 for all CPUs)")
    parser.add_argument("-v", "--value-types", nargs="+", default=None, help="Value types to analyze (e.g., espf espr). If not specified, all value types are analyzed.")
    parser.add_argument("--agg-levels", nargs="+", default=None, help="Aggregate levels to analyze: global, sequence, feature. If not specified, all levels are analyzed.")
    parser.add_argument("--feature-types", nargs="+", default=None, help="Feature types to analyze (e.g., gene exon RNA). If not specified, all feature types are analyzed.")
    parser.add_argument("--stat-test", default="nonparametric", 
                        choices=["auto", "parametric", "nonparametric", "welch", "kruskal"],
                        help="Statistical test selection: auto (test assumptions), parametric (Student/ANOVA), nonparametric (Mann-Whitney/Kruskal-Wallis, default), welch (Welch t-test/ANOVA)")
    parser.add_argument("--bmk-filter", nargs="+", default=["primary_padj", "kruskal_padj", "welch_padj", "anova_padj"],
                        help="Priority list of column names for filtering significant BMKs (default: primary_padj kruskal_padj welch_padj anova_padj). BMKs are collected in cascade until --max-bmks is reached.")
    parser.add_argument("--max-bmks", type=int, default=500,
                        help="Maximum number of BMKs to use for RandomForest and classification (default: 500). Prevents memory crashes with large datasets.")
    args = parser.parse_args()
    
    # Determine number of workers
    if args.jobs == -1:
        n_jobs = os.cpu_count()
    elif args.jobs < 1:
        parser.error("--jobs must be >= 1 or -1 for all CPUs")
    else:
        n_jobs = args.jobs
    
    log.info(f"Using {n_jobs} parallel job(s)")
    log.info(f"Statistical test method: {args.stat_test}")
    log.info(f"BMK filtering cascade: {' → '.join(args.bmk_filter)}")
    log.info(f"Max BMKs for ML models: {args.max_bmks}")

    # Check that at least one input file is provided
    if not args.aggregates and not args.features:
        parser.error("At least one of --aggregates or --features must be provided")

    outdir = args.outdir
    safe_mkdir(outdir)

    # Load data
    log.info("Loading data...")
    agg_df = None
    feat_df = None
    
    if args.aggregates:
        if not os.path.exists(args.aggregates):
            log.error(f"Aggregates file not found: {args.aggregates}")
            sys.exit(1)
        log.info(f"Loading aggregates from {args.aggregates}...")
        agg_df = pd.read_csv(args.aggregates, sep="\t", dtype={"SeqID": str, "Start": str, "End": str, "Strand": str})
        agg_df.columns = agg_df.columns.str.strip()  # Remove leading/trailing whitespace from column names
        # Strip whitespace from string columns
        for col in agg_df.select_dtypes(include=['object', 'string']).columns:
            agg_df[col] = agg_df[col].str.strip() if agg_df[col].dtype in ['object', 'string'] else agg_df[col]
        log.info(f"  Aggregates: {len(agg_df)} rows")
    else:
        log.info("Skipping aggregates (no file provided)")
    
    if args.features:
        if not os.path.exists(args.features):
            log.error(f"Features file not found: {args.features}")
            sys.exit(1)
        log.info(f"Loading features from {args.features}...")
        feat_df = pd.read_csv(args.features, sep="\t", dtype={"SeqID": str, "Start": str, "End": str, "Strand": str})
        feat_df.columns = feat_df.columns.str.strip()  # Remove leading/trailing whitespace from column names
        # Strip whitespace from string columns
        for col in feat_df.select_dtypes(include=['object', 'string']).columns:
            feat_df[col] = feat_df[col].str.strip() if feat_df[col].dtype in ['object', 'string'] else feat_df[col]
        log.info(f"  Features: {len(feat_df)} rows")
    else:
        log.info("Skipping features (no file provided)")

    # Parse sample information from whichever file is available
    sample_df = agg_df if agg_df is not None else feat_df
    sample_info = parse_sample_columns(sample_df.columns)
    all_value_types = get_value_types(sample_info)
    
    # Filter value types if specified
    if args.value_types:
        value_types = [vt for vt in args.value_types if vt in all_value_types]
        # Warn about non-existent value types
        missing = [vt for vt in args.value_types if vt not in all_value_types]
        if missing:
            log.warning(f"Requested value types not found in data: {missing}")
        if not value_types:
            log.error(f"None of the requested value types {args.value_types} found in data. Available: {all_value_types}")
            sys.exit(1)
    else:
        value_types = all_value_types
    
    # Count unique samples (group, sample, rep combinations)
    unique_samples = len(set((s["group"], s["sample"], s["rep"]) for s in sample_info))
    log.info(f"Available value types: {all_value_types}")
    if args.value_types:
        log.info(f"Analyzing value types: {value_types}")
    else:
        log.info(f"Analyzing all value types: {value_types}")
    log.info(f"Samples: {unique_samples} unique samples across {len(all_value_types)} value types ({len(sample_info)} total columns)")

    all_results = {}

    for vtype in value_types:
        log.info(f"\n{'='*60}")
        log.info(f"VALUE TYPE: {vtype}")
        log.info(f"{'='*60}")

        vtype_dir = os.path.join(outdir, vtype)
        safe_mkdir(vtype_dir)
        vcols = cols_for_vtype(sample_info, vtype)
        v_sample_info = sample_info_for_vtype(sample_info, vtype)
        all_results[vtype] = {"aggregate": {}, "feature": {}}
        
        # Collect all analysis tasks for this value_type
        tasks = []

        # ===============================================================
        # AGGREGATES
        # ===============================================================
        if agg_df is not None:
            log.info(f"\n--- AGGREGATES for {vtype} ---")
            agg_data = agg_df[agg_df["Mtype"] == "aggregate"].copy()

            agg_dir = os.path.join(vtype_dir, "aggregate")

            # Filter aggregate levels if specified
            agg_levels = args.agg_levels if args.agg_levels else ["global", "sequence", "feature"]
            log.info(f"Analyzing aggregate levels: {agg_levels}")

            # --- 1. Global aggregates (Type == global) ---
            if "global" in agg_levels:
                glob_agg = agg_data[agg_data["Type"] == "global"]
            else:
                glob_agg = pd.DataFrame()  # Empty dataframe if global is not requested
            
            if not glob_agg.empty:
                # All BMKs together (no Ptype/Ctype/Mode filter)
                tasks.append((
                    prepare_df_for_task(glob_agg), vcols, v_sample_info,
                    os.path.join(agg_dir, "global", "all_global_bmks"),
                    "Global - All BMKs",
                    ("aggregate", "global_all_bmks"),
                    args.stat_test,
                    args.bmk_filter,
                    args.max_bmks
                ))

                # all sites (Mode == all_sites)
                section = glob_agg[glob_agg["Mode"] == "all_sites"]
                tasks.append((
                    prepare_df_for_task(section), vcols, v_sample_info,
                    os.path.join(agg_dir, "global", "all_sites"),
                    "Global - All Sites",
                    ("aggregate", "global_all_sites"),
                    args.stat_test,
                    args.bmk_filter,
                    args.max_bmks
                ))

                # All sites by Ctype (Ptype == "." and 3 modes)
                for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                    section = glob_agg[(glob_agg["Ptype"] == ".") & (glob_agg["Mode"] == mode)]
                    tasks.append((
                        prepare_df_for_task(section), vcols, v_sample_info,
                        os.path.join(agg_dir, "global", f"by_ctype_{mode}"),
                        f"Global - By Ctype - {mode}",
                        ("aggregate", f"global_ctype_{mode}"),
                        args.stat_test,
                        args.bmk_filter,
                        args.max_bmks
                    ))

                # All sites by Ptype (Ptype != "." and 3 modes)
                ptypes = [p for p in glob_agg["Ptype"].unique() if p != "."]
                for ptype in ptypes:
                    for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                        section = glob_agg[(glob_agg["Ptype"] == ptype) & (glob_agg["Mode"] == mode)]
                        tasks.append((
                            prepare_df_for_task(section), vcols, v_sample_info,
                            os.path.join(agg_dir, "global", f"by_ptype_{ptype}_{mode}"),
                            f"Global - Ptype={ptype} - {mode}",
                            ("aggregate", f"global_ptype_{ptype}_{mode}"),
                            args.stat_test,
                            args.bmk_filter,
                            args.max_bmks
                        ))

            # --- 2. Chromosome/Sequence aggregates (Type == sequence) ---
            if "sequence" in agg_levels:
                seq_agg = agg_data[agg_data["Type"] == "sequence"]
            else:
                seq_agg = pd.DataFrame()  # Empty dataframe if chr is not requested
            
            if not seq_agg.empty:
                chromosomes = seq_agg["SeqID"].unique()
                for chrom in chromosomes:
                    chr_data = seq_agg[seq_agg["SeqID"] == chrom]

                    # All sites
                    section = chr_data[chr_data["Mode"] == "all_sites"]
                    tasks.append((
                        prepare_df_for_task(section), vcols, v_sample_info,
                        os.path.join(agg_dir, f"sequence/{chrom}", "all_sites"),
                        f"Chr {chrom} - All Sites",
                        ("aggregate", f"chr{chrom}_all_sites"),
                        args.stat_test,
                        args.bmk_filter,
                        args.max_bmks
                    ))

                    # By Ctype
                    for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                        section = chr_data[(chr_data["Ptype"] == ".") & (chr_data["Mode"] == mode)]
                        tasks.append((
                            prepare_df_for_task(section), vcols, v_sample_info,
                            os.path.join(agg_dir, f"sequence/{chrom}", f"by_ctype_{mode}"),
                            f"Chr {chrom} - By Ctype - {mode}",
                            ("aggregate", f"chr{chrom}_ctype_{mode}"),
                            args.stat_test,
                            args.bmk_filter,
                            args.max_bmks
                        ))

                    # By Ptype
                    local_ptypes = [p for p in chr_data["Ptype"].unique() if p != "."]
                    for ptype in local_ptypes:
                        for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                            section = chr_data[(chr_data["Ptype"] == ptype) & (chr_data["Mode"] == mode)]
                            tasks.append((
                                prepare_df_for_task(section), vcols, v_sample_info,
                                os.path.join(agg_dir, f"sequence/{chrom}", f"by_ptype_{ptype}_{mode}"),
                                f"Chr {chrom} - Ptype={ptype} - {mode}",
                                ("aggregate", f"chr{chrom}_ptype_{ptype}_{mode}"),
                                args.stat_test,
                                args.bmk_filter,
                                args.max_bmks
                            ))

                # --- All sequences combined (all chromosomes pooled) ---
                # all_sites
                section = seq_agg[seq_agg["Mode"] == "all_sites"]
                tasks.append((
                    prepare_df_for_task(section), vcols, v_sample_info,
                    os.path.join(agg_dir, "sequence", "all_sequence_bmks", "all_sites"),
                    "All Sequences - All Sites",
                    ("aggregate", "allseq_all_sites"),
                    args.stat_test,
                    args.bmk_filter,
                    args.max_bmks
                ))
                # by_ctype
                for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                    section = seq_agg[(seq_agg["Ptype"] == ".") & (seq_agg["Mode"] == mode)]
                    tasks.append((
                        prepare_df_for_task(section), vcols, v_sample_info,
                        os.path.join(agg_dir, "sequence", "all_sequence_bmks", f"by_ctype_{mode}"),
                        f"All Sequences - By Ctype - {mode}",
                        ("aggregate", f"allseq_ctype_{mode}"),
                        args.stat_test,
                        args.bmk_filter,
                        args.max_bmks
                    ))
                # by_ptype
                all_seq_ptypes = [p for p in seq_agg["Ptype"].unique() if p != "."]
                for ptype in all_seq_ptypes:
                    for mode in ["all_isoforms", "chimaera", "longest_isoform"]:
                        section = seq_agg[(seq_agg["Ptype"] == ptype) & (seq_agg["Mode"] == mode)]
                        tasks.append((
                            prepare_df_for_task(section), vcols, v_sample_info,
                            os.path.join(agg_dir, "sequence", "all_sequence_bmks", f"by_ptype_{ptype}_{mode}"),
                            f"All Sequences - Ptype={ptype} - {mode}",
                            ("aggregate", f"allseq_ptype_{ptype}_{mode}"),
                            args.stat_test,
                            args.bmk_filter,
                            args.max_bmks
                        ))

            # --- 3. Feature aggregates (Type == feature) ---
            if "feature" in agg_levels:
                feat_agg = agg_data[agg_data["Type"] == "feature"]
            else:
                feat_agg = pd.DataFrame()  # Empty dataframe if feature is not requested
            
            if not feat_agg.empty:
                fa_ptypes = feat_agg["Ptype"].unique()
                fa_ctypes = feat_agg["Ctype"].unique()
                fa_modes  = feat_agg["Mode"].unique()

                # --- all_feature_together: all chromosomes pooled ---
                for ptype in fa_ptypes:
                    for ctype in fa_ctypes:
                        for mode in fa_modes:
                            section = feat_agg[
                                (feat_agg["Ptype"] == ptype) &
                                (feat_agg["Ctype"] == ctype) &
                                (feat_agg["Mode"] == mode)
                            ]
                            safe_name = f"{ptype}_{ctype}_{mode}".replace(".", "all").replace("-", "_")
                            tasks.append((
                                prepare_df_for_task(section), vcols, v_sample_info,
                                os.path.join(agg_dir, "feature", "all_feature_together", safe_name),
                                f"Feature Agg - Ptype={ptype}, Ctype={ctype}, Mode={mode}",
                                ("aggregate", f"featagg_{safe_name}"),
                                args.stat_test,
                                args.bmk_filter,
                                args.max_bmks
                            ))

                # --- by_sequence: features grouped by chromosome ---
                fa_chroms = [c for c in feat_agg["SeqID"].unique() if c != "."]
                for chrom in fa_chroms:
                    chr_feat = feat_agg[feat_agg["SeqID"] == chrom]
                    chr_ptypes = chr_feat["Ptype"].unique()
                    chr_ctypes = chr_feat["Ctype"].unique()
                    chr_modes  = chr_feat["Mode"].unique()
                    for ptype in chr_ptypes:
                        for ctype in chr_ctypes:
                            for mode in chr_modes:
                                section = chr_feat[
                                    (chr_feat["Ptype"] == ptype) &
                                    (chr_feat["Ctype"] == ctype) &
                                    (chr_feat["Mode"] == mode)
                                ]
                                safe_name = f"{ptype}_{ctype}_{mode}".replace(".", "all").replace("-", "_")
                                tasks.append((
                                    prepare_df_for_task(section), vcols, v_sample_info,
                                    os.path.join(agg_dir, "feature", "by_sequence", str(chrom), safe_name),
                                    f"Feature Agg - Chr {chrom} - Ptype={ptype}, Ctype={ctype}, Mode={mode}",
                                    ("aggregate", f"featagg_chr{chrom}_{safe_name}"),
                                    args.stat_test,
                                    args.bmk_filter,
                                    args.max_bmks
                                ))
        else:
            log.info(f"\n--- Skipping AGGREGATES for {vtype} (no aggregates file) ---")

        # ===============================================================
        # FEATURES
        # ===============================================================
        if feat_df is not None:
            log.info(f"\n--- FEATURES for {vtype} ---")
            feat_data = feat_df.copy()  # Mtype is always "feature" in features file
            feat_dir = os.path.join(vtype_dir, "feature")

            # Build hierarchy
            tree, top_ids = build_feature_tree(feat_data)

            # Group top features by Type
            top_features = feat_data[feat_data["ParentIDs"] == "."]
            available_feature_types = top_features["Type"].unique().tolist()
            
            # Filter feature types if specified
            if args.feature_types:
                feature_types = [ft for ft in args.feature_types if ft in available_feature_types]
                missing = [ft for ft in args.feature_types if ft not in available_feature_types]
                if missing:
                    log.warning(f"Requested feature types not found in data: {missing}")
                if not feature_types:
                    log.error(f"None of the requested feature types {args.feature_types} found in data. Available: {available_feature_types}")
                    feature_types = []  # Will skip all features
            else:
                feature_types = available_feature_types
            
            if feature_types:
                log.info(f"Available feature types: {available_feature_types}")
                log.info(f"Analyzing feature types: {feature_types}")

            for ttype in feature_types:
                type_dir = os.path.join(feat_dir, ttype.replace(" ", "_"))
                type_features = top_features[top_features["Type"] == ttype]

                # Analyze all features of this type together
                all_ids_of_type = []
                for _, row in type_features.iterrows():
                    fid = row["ID"]
                    # Collect this feature and all descendants
                    def collect_ids(node_id):
                        ids = [node_id]
                        if node_id in tree:
                            for child in tree[node_id]["children"]:
                                ids.extend(collect_ids(child))
                        return ids
                    all_ids_of_type.extend(collect_ids(fid))

                type_all_df = feat_data[feat_data["ID"].isin(all_ids_of_type)]
                tasks.append((
                    prepare_df_for_task(type_all_df), vcols, v_sample_info,
                    os.path.join(type_dir, "_all"),
                    f"Features - {ttype} (all)",
                    ("feature", f"type_{ttype}_all"),
                    args.stat_test,
                    args.bmk_filter,
                    args.max_bmks
                ))

                # Per top-feature analysis
                for _, row in type_features.iterrows():
                    fid = row["ID"]
                    def collect_ids(node_id):
                        ids = [node_id]
                        if node_id in tree:
                            for child in tree[node_id]["children"]:
                                ids.extend(collect_ids(child))
                        return ids
                    sub_ids = collect_ids(fid)
                    sub_df = feat_data[feat_data["ID"].isin(sub_ids)]
                    safe_fid = fid.replace(":", "_").replace("/", "_")
                    tasks.append((
                        prepare_df_for_task(sub_df), vcols, v_sample_info,
                        os.path.join(type_dir, safe_fid),
                        f"Feature: {fid}",
                        ("feature", f"feature_{safe_fid}"),
                        args.stat_test,
                        args.bmk_filter,
                        args.max_bmks
                    ))
        else:
            log.info(f"\n--- Skipping FEATURES for {vtype} (no features file) ---")
        
        # Execute tasks (parallel or sequential)
        log.info(f"Executing {len(tasks)} analysis tasks...")
        if n_jobs > 1 and len(tasks) > 1:
            # Parallel execution with worker recycling to prevent memory leaks
            log.info(f"Submitting {len(tasks)} tasks to {n_jobs} workers...")
            
            # MEMORY FIX #1: Use spawn context to prevent fork Copy-on-Write memory inheritance
            # With fork (default on Linux), each worker inherits ALL parent memory via CoW
            # → 6 workers × 5GB parent = 30GB total (even if psutil shows less per process)
            # With spawn, workers start clean and only receive their task data via pickle
            mp_context = multiprocessing.get_context('spawn')
            log.info(f"  Using 'spawn' context (not fork) to avoid CoW memory inheritance")
            
            # Use max_tasks_per_child=20 to balance performance and memory
            # Higher than fork (was 5) because spawn workers are clean but have import overhead (~275MB)
            log.info(f"  Worker recycling: Each worker will process max 20 tasks before restart")
            
            # CRITICAL: Use lazy submission to avoid loading all 1491 tasks (dict copies) in memory at once
            # Submit only max_pending_tasks at a time, then submit new ones as they complete
            max_pending_tasks = n_jobs * 3  # Keep 3x workers worth of tasks in flight
            log.info(f"  Lazy submission: Max {max_pending_tasks} tasks in memory at once (was {len(tasks)})")
            
            with ProcessPoolExecutor(max_workers=n_jobs, max_tasks_per_child=20, mp_context=mp_context) as executor:
                future_to_key = {}
                task_iter = iter(tasks)
                completed = 0
                failed = 0
                submitted_count = 0
                last_progress_time = time.time()
                stall_timeout = 90  # Increased from 40s for spawn startup delay (imports take ~2-5s per worker)
                last_worker_check = time.time()
                worker_check_interval = 15  # Check worker health every 15 seconds
                
                # Initial submission of first batch
                initial_batch = min(max_pending_tasks, len(tasks))
                log.info(f"  Submitting initial batch of {initial_batch} tasks...")
                for _ in range(initial_batch):
                    try:
                        df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks = next(task_iter)
                        future = executor.submit(analyze_section_wrapper, (df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks))
                        future_to_key[future] = (key, name)
                        submitted_count += 1
                        # Explicitly free the df_source reference to help GC
                        del df_source
                    except StopIteration:
                        break
                
                log.info(f"  Initial batch submitted. Processing with {n_jobs} workers...")
                log.info(f"  Note: Tasks have a 2-minute timeout. Anti-deadlock active.")
                
                # Process futures with stall detection
                pending_futures = set(future_to_key.keys())
                
                while pending_futures:
                    # Periodically check if workers have died
                    current_time = time.time()
                    if psutil is not None and (current_time - last_worker_check) > worker_check_interval:
                        try:
                            process = psutil.Process()
                            children = process.children(recursive=True)
                            n_active_workers = len(children)
                            
                            # If more than half of workers are dead, we have a problem
                            if n_active_workers < (n_jobs / 2) and len(pending_futures) > n_active_workers * 2:
                                log.error(f"  ✗ WORKER DEATH DETECTED: Only {n_active_workers}/{n_jobs} workers alive with {len(pending_futures)} pending tasks")
                                log.error(f"     Most likely cause: Out-Of-Memory (OOM) killed workers")
                                log.error(f"     Cancelling all pending tasks to prevent infinite hang")
                                log.error(f"     TIP: Restart with fewer workers (--n-jobs 2 or --n-jobs 3)")
                                
                                # Cancel all pending futures
                                for future in list(pending_futures):
                                    future.cancel()
                                    key, section_name = future_to_key[future]
                                    log.error(f"     ✗ CANCELLED (orphaned): {section_name}")
                                    failed += 1
                                break
                            
                            last_worker_check = current_time
                        except Exception as e:
                            log.debug(f"Worker health check failed: {e}")
                    
                    # Use short timeout on as_completed to check for stalls
                    try:
                        done_iter = as_completed(pending_futures, timeout=5)
                        for future in done_iter:
                            pending_futures.discard(future)
                            key, section_name = future_to_key[future]
                            
                            # Process the completed future
                            try:
                                result_key, results = future.result(timeout=1)  # Short timeout since already done
                                mtype, section_key = result_key
                                all_results[vtype][mtype][section_key] = results
                                completed += 1
                                last_progress_time = time.time()  # Reset stall timer
                                
                                # MEMORY FIX #3: Immediate cleanup after storing (results is slim dict with only paths)
                                del results
                                
                                # Force garbage collection every 10 tasks to free memory in main process
                                if completed % 10 == 0:
                                    gc.collect()
                                
                                # Progress update every 5 completions or at key milestones
                                if completed % 5 == 0 or completed in [1, 10, 25, 50, 100]:
                                    log.info(f"  ✓ Progress: {completed}/{len(tasks)} completed, {failed} failed, {submitted_count - completed - failed} submitted pending")
                            except TimeoutError:
                                failed += 1
                                last_progress_time = time.time()
                                log.error(f"  ✗ TIMEOUT ({completed+failed}/{len(tasks)}): {section_name} - exceeded 2 minutes")
                            except Exception as e:
                                failed += 1
                                last_progress_time = time.time()
                                # Check if it's a worker crash (common patterns in error message)
                                if "process" in str(e).lower() and ("terminate" in str(e).lower() or "crash" in str(e).lower() or "abrupt" in str(e).lower()):
                                    log.error(f"  ✗ CRASH ({completed+failed}/{len(tasks)}): {section_name} - worker OOM or crash")
                                else:
                                    log.error(f"  ✗ ERROR ({completed+failed}/{len(tasks)}): {section_name} - {type(e).__name__}")
                            
                            # CRITICAL: Submit next task to maintain max_pending_tasks in flight
                            # This lazy submission keeps memory usage constant regardless of total tasks
                            if len(pending_futures) < max_pending_tasks:
                                try:
                                    df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks = next(task_iter)
                                    new_future = executor.submit(analyze_section_wrapper, (df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks))
                                    future_to_key[new_future] = (key, name)
                                    pending_futures.add(new_future)
                                    submitted_count += 1
                                    # Explicitly free the df_source reference to help GC
                                    del df_source
                                    # Log every 50 submissions
                                    if submitted_count % 50 == 0:
                                        log.info(f"  → Submitted up to {submitted_count}/{len(tasks)} tasks (lazy mode)")
                                except StopIteration:
                                    pass  # No more tasks to submit
                            
                            # Break inner loop to check stall timeout
                            break
                    except TimeoutError:
                        # No futures completed in 5 seconds, check for stall
                        elapsed_since_progress = time.time() - last_progress_time
                        if elapsed_since_progress > stall_timeout:
                            log.error(f"  ✗ DEADLOCK DETECTED: No progress for {elapsed_since_progress:.0f}s")
                            log.error(f"     Cancelling {len(pending_futures)} remaining tasks to prevent infinite hang")
                            # Cancel all pending futures
                            for future in pending_futures:
                                future.cancel()
                                key, section_name = future_to_key[future]
                                log.error(f"     ✗ CANCELLED: {section_name}")
                                failed += 1
                            break
                
                if failed > 0:
                    log.warning(f"  {failed} tasks failed or cancelled out of {len(tasks)} total")
        else:
            # Sequential execution
            for df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks in tasks:
                # Reconstruct DataFrame from pickled format
                if isinstance(df_source, tuple) and df_source[0] == 'pickled':
                    import pickle
                    df = pickle.loads(df_source[1])
                else:
                    df = pd.DataFrame(df_source)
                results = analyze_section(df, cols, info, outdir, name, stat_test=stat_test, bmk_filter_cols=bmk_filter, max_bmks=max_bmks)
                
                # Create slim results dict (same as parallel mode for consistency)
                slim_results = {
                    "differential_table": os.path.join(outdir, "differential", "differential_results.csv")
                }
                
                mtype, section_key = key
                all_results[vtype][mtype][section_key] = slim_results

    # ===============================================================
    # GLOBAL RANKING
    # ===============================================================
    log.info("\n--- GLOBAL RANKING ---")
    global_ranking(all_results, os.path.join(outdir, "global_ranking"))

    # Save manifest
    manifest = {
        "value_types": value_types,
        "outdir": outdir,
        "n_aggregates": len(agg_df) if agg_df is not None else 0,
        "n_features": len(feat_df) if feat_df is not None else 0,
        "sample_info": sample_info,
    }
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    log.info(f"\nAnalysis complete. Results saved to {outdir}/")


if __name__ == "__main__":
    main()
