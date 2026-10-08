"""barometer.analysis – Per-section statistical / ML analyses.

Each function analyzes one aspect of a section (QC, batch effect, descriptive
stats, multivariate, differential, correlation, ranking, classification,
stability, heatmap) or orchestrates them (analyze_section, build_feature_tree).
"""

import logging
import os
import subprocess
import tempfile
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import StandardScaler, LabelEncoder
from statsmodels.stats.multitest import multipletests

from .utils import (
    safe_mkdir,
    numeric_df,
    save_fig,
    group_for_col,
)

log = logging.getLogger(__name__)


def append_res_df_by_df_via_index(res_df, df, sample_cols, meta_cols_exclude):
    """Join metadata columns from df onto res_df (aligned by index), metadata first."""
    meta_cols = [
        c for c in df.columns
        if c not in sample_cols
        and c not in meta_cols_exclude
        and not c.endswith(("::successes", "::trials"))
    ]
    # Merge with original df to get metadata
    res_df = res_df.join(df.loc[res_df.index, meta_cols], how="left")
    # ⭐ RÉORGANISER : mettre les METADATA en PREMIER ⭐
    existing_meta = [c for c in meta_cols if c in res_df.columns]
    stat_cols = [c for c in res_df.columns if c not in meta_cols]
    # Ordre final : metadata d'abord, puis statistiques
    res_df = res_df[existing_meta + stat_cols]
    return res_df


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


def run_beta_binomial_analysis(df, sample_info, groups, outdir):
    """Fit beta-binomial models using DRIP's unrounded success/trial columns."""
    pair_keys = [(g1, g2, f"{g1}_vs_{g2}") for g1, g2 in combinations(groups, 2)]
    output = pd.DataFrame(index=df.index)
    output["beta_binomial_status"] = "no_covered_samples"
    output["beta_binomial_stat"] = np.nan
    output["beta_binomial_pval"] = np.nan
    output["beta_binomial_dispersion"] = np.nan
    for group in groups:
        output[f"beta_binomial_mean_{group}"] = np.nan
    for _, _, pair_key in pair_keys:
        output[f"beta_binomial_log_odds_{pair_key}"] = np.nan
        output[f"beta_binomial_diff_{pair_key}"] = np.nan
        output[f"beta_binomial_pval_{pair_key}"] = np.nan

    count_frames = []
    for sample in sample_info:
        success_col = f"{sample['col']}::successes"
        trials_col = f"{sample['col']}::trials"
        if success_col not in df.columns or trials_col not in df.columns:
            raise ValueError(
                f"Missing beta-binomial counts for {sample['col']}; "
                "rerun DRIP with count columns enabled."
            )

        successes = pd.to_numeric(df[success_col], errors="coerce")
        trials = pd.to_numeric(df[trials_col], errors="coerce")
        valid = (
            successes.notna()
            & trials.notna()
            & (trials > 0)
            & (successes >= 0)
            & (successes <= trials)
        )
        if valid.any():
            count_frames.append(pd.DataFrame({
                "row_id": df.index[valid].map(str),
                "group": sample["group"],
                "sample": sample["sample"],
                "replicate": sample["rep"],
                "successes": successes.loc[valid].astype(int).values,
                "trials": trials.loc[valid].astype(int).values,
            }))

    if not count_frames:
        return output

    counts = pd.concat(count_frames, ignore_index=True)
    # The R script lives alongside this module in bin/barometer/.
    # The env var still takes priority.
    script_path = Path(os.environ.get(
        "RAIN_BETA_BINOMIAL_SCRIPT",
        Path(__file__).resolve().with_name("barometer_beta_binomial.R"),
    ))
    if not script_path.exists():
        raise FileNotFoundError(f"Beta-binomial R script not found: {script_path}")

    with tempfile.TemporaryDirectory(prefix="barometer_beta_binomial_") as temp_dir:
        input_path = Path(temp_dir) / "counts.tsv"
        output_path = Path(temp_dir) / "results.tsv"
        counts.to_csv(input_path, sep="\t", index=False)
        try:
            completed = subprocess.run(
                ["Rscript", str(script_path), str(input_path), str(output_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                "Beta-binomial model failed: "
                f"{exc.stderr[-4000:]}"
            ) from exc

        if not output_path.exists():
            raise RuntimeError(
                "Beta-binomial script did not create a results file. "
                f"{completed.stderr[-4000:]}"
            )
        fitted = pd.read_csv(output_path, sep="\t", dtype={"row_id": str})

    index_by_row_id = {str(index): index for index in df.index}
    for row in fitted.itertuples(index=False):
        index = index_by_row_id.get(str(row.row_id))
        if index is None:
            continue
        if row.test == "global":
            output.at[index, "beta_binomial_stat"] = row.statistic
            output.at[index, "beta_binomial_pval"] = row.p_value
            output.at[index, "beta_binomial_dispersion"] = row.dispersion
            output.at[index, "beta_binomial_status"] = row.status
        elif row.test == "mean":
            output.at[index, f"beta_binomial_mean_{row.group1}"] = row.mean1
        elif row.test == "pairwise":
            pair_key = f"{row.group1}_vs_{row.group2}"
            if pair_key in {key for _, _, key in pair_keys}:
                output.at[index, f"beta_binomial_log_odds_{pair_key}"] = row.estimate
                output.at[index, f"beta_binomial_diff_{pair_key}"] = row.mean1 - row.mean2
                output.at[index, f"beta_binomial_pval_{pair_key}"] = row.p_value

    return output


def differential_analysis(df, sample_cols, sample_info, outdir, stat_test="auto"):
    """Pairwise and global differential tests between groups.

    Args:
        stat_test: 'auto', 'parametric', 'nonparametric', 'welch', 'kruskal', 'beta-binomial'
            - auto: test normality/homogeneity and choose best test
            - parametric: Student t-test / ANOVA (assumes normality + equal variances)
            - nonparametric: Mann-Whitney U / Kruskal-Wallis (no assumptions)
            - welch: Welch t-test / Welch ANOVA (assumes normality, unequal variances OK)
            - kruskal: alias for nonparametric (backward compatibility)
            - beta-binomial: model the DRIP success/trial COUNTS directly with a
              beta-binomial GLMM (see below).

    Beta-binomial mode (stat_test == "beta-binomial")
    -------------------------------------------------
    The other modes test the *rounded proportions* (the aggregate values in the
    sample columns). The beta-binomial mode instead uses the raw, unrounded
    counts that DRIP emits per sample:

        {sample}::successes  -> number of reads showing the edited base
        {sample}::trials     -> total number of reads covering the site

    For each biomarker a single beta-binomial GLMM is fitted (via R/glmmTMB,
    script selected by the RAIN_BETA_BINOMIAL_SCRIPT env var, default
    barometer_beta_binomial.R):

        cbind(successes, trials - successes) ~ condition

    The beta-binomial family adds a dispersion parameter on top of the
    binomial, which absorbs the extra-binomial overdispersion that is typical
    of sequencing/editing data (reads are not independent Bernoulli trials).
    This makes the p-values far more reliable than a plain binomial or a test
    on rounded proportions.

    From that single model two kinds of Wald tests are derived:

      * GLOBAL test  (1 per BMK): "does condition have ANY effect?"
          -> beta_binomial_stat / beta_binomial_pval
             -> beta_binomial_padj (FDR BH)
          This is the beta-binomial equivalent of an ANOVA F-test.

      * PAIRWISE tests (C(n_groups, 2) per BMK): one Wald contrast per pair
          -> beta_binomial_pval_{g1}_vs_{g2}
             -> beta_binomial_padj_{g1}_vs_{g2} (FDR BH)

    A BMK can be significant in one pairwise contrast while the global test is
    not (the effect is concentrated on a single pair and diluted over the
    global degrees of freedom), and vice-versa. Significance is therefore
    defined as "padj < 0.05 in ANY test (global or pairwise)", consistent with
    the section-level logic.

    In this mode the primary_* columns are populated from the beta-binomial
    results (primary_padj = beta_binomial_padj, primary_padj_{pair} =
    beta_binomial_padj_{pair}), and the near-zero-variance prefilter is
    skipped (counts, not proportions, drive the inference).

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
            # DEBUG: Log first BMK's pairwise test setup
            if idx == ndf.index[0]:
                pairs_list = list(combinations(groups, 2))
                log.debug(f"PAIRWISE: {len(groups)} groups → {len(pairs_list)} pairs")
                log.debug(f"PAIRWISE: Groups = {groups}")
                log.debug(f"PAIRWISE: Pairs = {pairs_list}")
            
            for g1, g2 in combinations(groups, 2):
                v1, v2 = group_values.get(g1, []), group_values.get(g2, [])
                pair_key = f"{g1}_vs_{g2}"
                
                # DEBUG: Log first BMK's pairwise data
                if idx == ndf.index[0]:
                    log.debug(f"PAIRWISE: {pair_key} → v1={len(v1)} v2={len(v2)}")
                
                if len(v1) >= 1 and len(v2) >= 1:
                    # DEBUG: Log first BMK entering pairwise test block
                    if idx == ndf.index[0]:
                        log.debug(f"PAIRWISE: Entering test block for {pair_key}")
                    
                    # Mann-Whitney U (nonparametric, always compute)
                    try:
                        stat, pval = stats.mannwhitneyu(v1, v2, alternative="two-sided")
                        row_data[f"mwu_stat_{pair_key}"] = stat
                        row_data[f"mwu_pval_{pair_key}"] = pval
                        
                        # DEBUG: Confirm keys added
                        if idx == ndf.index[0]:
                            log.debug(f"PAIRWISE: Added mwu keys for {pair_key}: stat={stat:.3f}, pval={pval:.3e}")
                    except Exception as e:
                        row_data[f"mwu_stat_{pair_key}"] = np.nan
                        row_data[f"mwu_pval_{pair_key}"] = np.nan
                        if idx == ndf.index[0]:
                            log.debug(f"PAIRWISE: Mann-Whitney failed for {pair_key}: {e}")
                    
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

                # Ajout des colonnes primary_stat/pval/padj pour chaque comparaison
                # On détermine le test principal pour cette paire
                if test_choice == "parametric":
                    stat = row_data.get(f"student_stat_{pair_key}", np.nan)
                    pval = row_data.get(f"student_pval_{pair_key}", np.nan)
                    padj = row_data.get(f"student_padj_{pair_key}", np.nan)
                    test_used = "student"
                elif test_choice == "welch":
                    stat = row_data.get(f"welch_stat_{pair_key}", np.nan)
                    pval = row_data.get(f"welch_pval_{pair_key}", np.nan)
                    padj = row_data.get(f"welch_padj_{pair_key}", np.nan)
                    test_used = "welch"
                else:  # nonparametric
                    stat = row_data.get(f"mwu_stat_{pair_key}", np.nan)
                    pval = row_data.get(f"mwu_pval_{pair_key}", np.nan)
                    padj = row_data.get(f"mwu_padj_{pair_key}", np.nan)
                    test_used = "mwu"
                row_data[f"primary_stat_{pair_key}"] = stat
                row_data[f"primary_pval_{pair_key}"] = pval
                row_data[f"primary_test_{pair_key}"] = test_used
    
            rows.append(row_data)

    # Convert rows to DataFrame
    # DEBUG: Check last row too to ensure all rows have pairwise columns
    if len(rows) > 1:
        last_row_keys = list(rows[-1].keys())
        pval_keys_last = [k for k in last_row_keys if k.endswith("_pval")]
        log.debug(f"Last row has {len(last_row_keys)} total keys")
        log.debug(f"Last row has {len(pval_keys_last)} _pval keys")
    
    # Convert rows to DataFrame
    res_df = pd.DataFrame(rows)
    res_df = res_df.set_index("index")

    if stat_test == "beta-binomial":
        beta_results = run_beta_binomial_analysis(df, sample_info, groups, outdir)
        res_df = res_df.join(beta_results)
    
    # ============================================================================
    # FDR CORRECTION - Multiple testing correction (FDR Benjamini-Hochberg) for ALL p-value columns
    # ============================================================================
    log.info(f"  Applying FDR correction to p-values...")

    pval_cols = [c for c in res_df.columns if "pval" in str(c).lower()]

    # Adding the padj columns in-place may emit a harmless pandas
    # PerformanceWarning ("DataFrame is fragmented")
    n_padj_created = 0
    for col in pval_cols:  # FIXED: Iterate over filtered list instead of checking endswith
        adj_col = col.replace("_pval", "_padj")
        res_df[adj_col] = np.nan  # Always create the column

        pvals = res_df[col].values
        mask = ~np.isnan(pvals)
        if mask.sum() > 0:
            _, corrected, _, _ = multipletests(pvals[mask], method="fdr_bh")
            res_df.loc[mask, adj_col] = corrected
        n_padj_created += 1
    log.info(f"  Created {n_padj_created} adjusted p-value columns")

    if stat_test == "beta-binomial":
        res_df["primary_stat"] = res_df["beta_binomial_stat"]
        res_df["primary_pval"] = res_df["beta_binomial_pval"]
        res_df["primary_padj"] = res_df["beta_binomial_padj"]
        for g1, g2 in combinations(groups, 2):
            pair_key = f"{g1}_vs_{g2}"
            res_df[f"diff_{pair_key}"] = res_df[f"beta_binomial_diff_{pair_key}"]
            res_df[f"primary_stat_{pair_key}"] = res_df[f"beta_binomial_log_odds_{pair_key}"]
            res_df[f"primary_pval_{pair_key}"] = res_df[f"beta_binomial_pval_{pair_key}"]
            res_df[f"primary_padj_{pair_key}"] = res_df[f"beta_binomial_padj_{pair_key}"]
            res_df[f"primary_test_{pair_key}"] = "beta-binomial"

    # ============================================================================
    # ADD ORIGINAL METADATA from input df
    # ============================================================================
    meta_cols_exclude = ["some_col", "another_col"]  # Example of columns to exclude from metadata (if any)
    res_df = append_res_df_by_df_via_index(res_df, df, sample_cols, meta_cols_exclude)  # Ensure index is aligned for joining

    # ============================================================================
    # OUTPUT RESULTS
    # ============================================================================
    res_df.to_csv(os.path.join(outdir, "differential_results.csv"), index=True)
    results["table"] = os.path.join(outdir, "differential_results.csv")
    results["stat_test_method"] = stat_test

    # Extract significant biomarkers (padj < 0.05 in ANY test)
    if stat_test == "beta-binomial":
        padj_cols = [c for c in res_df.columns
                     if c == "beta_binomial_padj" or c.startswith("beta_binomial_padj_")]
    else:
        padj_cols = [c for c in res_df.columns if c.endswith("_padj")]
    if padj_cols:
        # Create boolean mask: True if ANY padj column is < 0.05
        sig_mask = res_df[padj_cols].lt(0.05).any(axis=1)
        sig_df = res_df[sig_mask].copy()
        
        if len(sig_df) > 0:
            # Add summary column: count of significant tests per BMK
            sig_df["n_significant_tests"] = res_df[padj_cols].lt(0.05).sum(axis=1)[sig_mask]
            
            # Sort by minimum padj value (most significant first)
            sig_df["min_padj"] = res_df[padj_cols].min(axis=1)[sig_mask]
            sig_df = sig_df.sort_values("min_padj")
            
            sig_file = os.path.join(outdir, "significant_biomarkers.csv")
            sig_df.to_csv(sig_file, index=False)
            results["significant_table"] = sig_file
            log.info(f"  Found {len(sig_df)} significant biomarkers (padj < 0.05 in any test)")
        else:
            log.info(f"  No significant biomarkers found (padj < 0.05)")
    
    # Volcano-like plot for each pairwise comparison (using primary test)
    log.info(f"  Generating volcano plots for {len(list(combinations(groups, 2)))} pairwise comparisons...")
    n_plots_created = 0
    for g1, g2 in combinations(groups, 2):
        pair_key = f"{g1}_vs_{g2}"
        diff_col = f"diff_{pair_key}"
        
        # Choose p-value column based on test method
        if stat_test == "beta-binomial":
            diff_col = f"beta_binomial_diff_{pair_key}"
            pval_col = f"beta_binomial_pval_{pair_key}"
            padj_col = f"beta_binomial_padj_{pair_key}"
            test_label = "Beta-binomial"
        elif stat_test == "parametric":
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
        
        # Check if columns exist
        if diff_col not in res_df.columns:
            log.warning(f"  Column {diff_col} not found, skipping volcano plot for {pair_key}")
            continue
        if padj_col not in res_df.columns:
            log.warning(f"  Column {padj_col} not found, skipping volcano plot for {pair_key}")
            continue
            
        pdf = res_df[[diff_col, padj_col]].dropna()
        if len(pdf) > 0:
            fig, ax = plt.subplots(figsize=(8, 6))
            neg_log_p = -np.log10(pdf[padj_col].clip(lower=1e-300))
            
            # Color points by significance level
            colors = []
            for p in pdf[padj_col]:
                if p < 0.05:
                    colors.append("red")
                elif p < 0.1:
                    colors.append("orange")
                else:
                    colors.append("lightgrey")
            
            ax.scatter(pdf[diff_col], neg_log_p, c=colors, alpha=0.6, s=20, edgecolors='none')
            
            # Add significance threshold lines
            ax.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=1.5, alpha=0.7, label="p = 0.05")
            ax.axhline(-np.log10(0.1), color="green", linestyle="--", linewidth=1.5, alpha=0.7, label="p = 0.1")
            
            # Add legend for significance levels
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='red', alpha=0.6, label='p < 0.05'),
                Patch(facecolor='orange', alpha=0.6, label='0.05 ≤ p < 0.1'),
                Patch(facecolor='lightgrey', alpha=0.6, label='p ≥ 0.1')
            ]
            legend1 = ax.legend(handles=legend_elements, loc='upper left', title='Significativité')
            ax.add_artist(legend1)
            
            # Add legend for threshold lines
            ax.legend(loc='upper right', title='Seuils')
            
            ax.set_xlabel(f"Différence de moyenne ({g1} - {g2})")
            ax.set_ylabel("-log10(p-value ajustée)")
            ax.set_title(f"Volcano plot : {g1} vs {g2} ({test_label})")
            ax.grid(True, alpha=0.3, linestyle=':')
            save_fig(fig, os.path.join(outdir, f"volcano_{pair_key}.png"))
            n_plots_created += 1
            log.info(f"    Created volcano plot : {pair_key}")
        else:
            log.warning(f"  No data available for volcano plot {pair_key}")
    
    log.info(f"  Generated {n_plots_created} volcano plots")

    return results

# ---------------------------------------------------------------------------
# Correlation / Network Analysis - with clustering
# ---------------------------------------------------------------------------

def correlation_network(df, sample_cols, outdir, max_bmks=50):
    """Compute BMK-BMK correlation matrix and heatmap with hierarchical clustering.
    
    MEMORY FIX: Reduced max_bmks from 200 to 50 to prevent OOM.
    - Correlation matrix: N×N float64 = N² × 8 bytes
    - With N=200: 320 KB per section × 1491 sections = 477 MB accumulated
    - With N=50: 20 KB per section × 1491 sections = 30 MB accumulated (16x reduction)
    - Heatmap figure with N=50: ~4 MB instead of 64 MB per image
    """
    safe_mkdir(outdir)
    results = {}
    
    df = df.set_index("uid", drop=False)

    # avoid Correlation analysis failed: could not convert string to float: 'NA'
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat = mat.dropna(axis=0, how="all")
    
    # Early exit if too few BMKs
    if len(mat) < 3:
        log.info(f"  Skipping correlation analysis: only {len(mat)} BMKs (need at least 3)")
        return results
    
    var = mat.var(axis=1)
    top_idx = var.nlargest(min(max_bmks, len(var))).index
    sub = mat.loc[top_idx]

    if len(sub) >= 2:
        # Compute correlation matrix (memory intensive: N×N)
        corr = sub.T.corr()
        corr.to_csv(os.path.join(outdir, "bmk_correlation.csv"))

        # ── Heatmap simple (sans clustering) ──────────────────────────────────
        fig_width  = min(45, max(16, len(sub) * 0.40))  
        fig_height = min(40, max(12, len(sub) * 0.40))   

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(corr, cmap="RdBu_r", center=0, ax=ax,
                    xticklabels=len(sub) < 50,
                    yticklabels=len(sub) < 50)
        ax.set_title(f"BMK Correlation (top {len(sub)} by variance)")
        save_fig(fig, os.path.join(outdir, "bmk_correlation_heatmap.png"))

        # ── Clustermap (regroupement par ressemblance) ────────────────────────
        # sns.clustermap retourne une ClusterGrid, pas un Figure+Axes classique
        if len(sub) >= 3:  # clustering nécessite au moins 3 éléments
            try:
                cmap_width  = min(45, max(16,  len(sub) * 0.40))
                cmap_height = min(45, max(16,  len(sub) * 0.40))

                cg = sns.clustermap(
                    corr,
                    cmap="RdBu_r",
                    center=0,
                    method="average",       # liaison hiérarchique : average, complete, ward…
                    metric="euclidean",     # distance entre lignes/colonnes
                    figsize=(cmap_width, cmap_height),
                    xticklabels=len(sub) < 50,
                    yticklabels=len(sub) < 50,
                    dendrogram_ratio=0.15,  # taille relative du dendrogramme
                    cbar_pos=(0.02, 0.8, 0.03, 0.18),  # position colorbar
                )
                cg.fig.suptitle(
                    f"BMK Correlation – clustered (top {len(sub)} by variance)",
                    y=1.02
                )

                # Sauvegarder via la figure de la ClusterGrid
                save_fig(cg.fig, os.path.join(outdir, "bmk_correlation_clustermap.png"))

                # Ordre des BMKs après clustering (utile pour debug / export)
                clustered_order = corr.index[cg.dendrogram_row.reordered_ind].tolist()
                log.debug(f"Clustered BMK order: {clustered_order[:10]}{'...' if len(clustered_order) > 10 else ''}")

            except Exception as e:
                log.warning(f"  Clustermap failed (skipping): {e}")

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

    # ⭐ CRÉER bmk_ids BASÉ SUR ndf.index ⭐
    if "uid" in df.columns:
        bmk_ids = df.loc[ndf.index, "uid"].values
    else:
        bmk_ids = ndf.index.values

    # Build matrix: rows = samples, cols = BMKs
    mat = ndf[sample_cols].T.copy()

    # Use the unique index (already set in analyze_section)
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
    diff_file = os.path.join(os.path.dirname(outdir), "5_differential", "differential_results.csv")
    selected_bmks = set()
    n_bmks_total = len(mat.columns)
    filter_log = []
    
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        
        # ⭐ Déterminer quelle colonne utiliser pour l'ID ⭐
        id_col = "uid" if "uid" in diff_df.columns else "index" if "index" in diff_df.columns else None
        
        if id_col is None:
            log.warning("No 'uid' or 'index' column found in differential results")
        else:
            # PASS 1: Cascade through adjusted p-values (strong evidence)
            for filter_col in bmk_filter_cols:
                if len(selected_bmks) >= max_bmks:
                    break
                    
                if filter_col in diff_df.columns:
                    sig_bmks = diff_df[diff_df[filter_col] < 0.05][id_col].values
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
                    if pval_col != filter_col and pval_col in diff_df.columns :
                        sig_bmks = diff_df[diff_df[pval_col] < 0.05][id_col].values
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
    diff_file = os.path.join(os.path.dirname(outdir), "5_differential", "differential_results.csv")
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        id_col = "uid" if "uid" in diff_df.columns else "index" if "index" in diff_df.columns else None

        if "kruskal_padj" in diff_df.columns:
            rank_df = diff_df[[id_col, "kruskal_padj"]].dropna().sort_values("kruskal_padj")
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

    bmk_ids = df["uid"].values if "uid" in df.columns else [str(i) for i in range(len(df))]

    mat = ndf[sample_cols].T.fillna(0).copy()
    
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
    diff_file = os.path.join(os.path.dirname(outdir), "5_differential", "differential_results.csv")
    selected_bmks = set()
    filter_log = []
    
    if os.path.exists(diff_file):
        diff_df = pd.read_csv(diff_file)
        
        # Cascade through filter columns in priority order
        for filter_col in bmk_filter_cols:
            if len(selected_bmks) >= max_bmks:
                break  # Stop if we have enough
                
            if filter_col in diff_df.columns and "uid" in diff_df.columns:
                # Get significant BMKs from this column
                sig_bmks = diff_df[diff_df[filter_col] < 0.05]["uid"].values
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
    else:
        log.warning(f"Go ahead we have enough biomarkers for classification")

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
    cv_df = pd.DataFrame({"bmk": df["uid"].values if "uid" in df.columns else range(len(df)),
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
            elif "uid" in df.columns:
                base_label = str(df.loc[idx, "uid"])
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
    save_fig(fig, os.path.join(outdir, "10_heatmap.png"))


# ---------------------------------------------------------------------------
# Top-level orchestration per section
# ---------------------------------------------------------------------------
# df arrive complete here
# sample_cols looks like'condition::sample::replicate::valuetype'
def analyze_section(df, sample_cols, sample_info, outdir, section_name, stat_test="auto", bmk_filter_cols=None, max_bmks=500, enabled_tests=None):
    """Run all analyses on a filtered dataframe for a given report section.
    
    enabled_tests : list[str] | None
        - None ou liste vide  → tous les tests sont exécutés
        - ["qc", "pca", ...]  → seuls les tests listés sont exécutés
    
    Noms valides :
        "filter1", "filter2",
        "qc", "batch", "descriptive", "multivariate", "differential",
        "correlation", "ranking", "classification", "stability", "heatmap"
        
    info:
        filter1 = near-zero variance filter
        filter2 = size limit filter (top 10k)
    """

    # ── Résolution de la liste active ────────────────────────────────────────
    ALL_TESTS = ["filter1", "filter2","qc", "batch", "descriptive", "multivariate", "differential","correlation", "ranking", "classification", "stability", "heatmap"]
    
    # Check size for skipping empty sections
    len_df = len(df)
    len_sample_cols = len(sample_cols)
    if len_df == 0:
        log.info(f"  Section '{section_name}' is empty, skipping.")
        return {}

    # Deal with tests
    if not enabled_tests:          # None ou liste vide → tout activer
        active = set(ALL_TESTS)
    else:
        active = set(enabled_tests)
        unknown = active - set(ALL_TESTS)
        if unknown:
            log.warning(f"Unknown test names ignored: {unknown}")
        active &= set(ALL_TESTS)   # on ne garde que les noms valides

    log.info(f"Active analyses: {sorted(active)}")

    # ── Helper pour sauter proprement un test désactivé ──────────────────────
    def is_active(name):
        if name in active:
            return True
        log.info(f"Skipping '{name}' (not in enabled_tests)")
        return False


    # Create log prefix for this section
    log_prefix = f"[{section_name[:50]}]"  # Limit to 50 chars
    
    log.info(f"{log_prefix} Analyzing section: {section_name} ({len_df} BMKs, {len_sample_cols} samples)")
    safe_mkdir(outdir)    
    results = {"n_bmks": len_df, "n_samples": len_sample_cols, "section": section_name}
    
    if is_active("filter1") and stat_test != "beta-binomial":
        # PRE-FILTER #1: Remove BMKs with near-zero variance FIRST (before any size limit)
        # This prevents numerical issues and removes uninformative BMKs early
        ndf_prefilter = numeric_df(df, sample_cols)
        variance_prefilter = ndf_prefilter[sample_cols].var(axis=1)
        variable_mask = variance_prefilter > 1e-10
        n_constants = (~variable_mask).sum()
        
        if n_constants > 0:
            log.info(f"{log_prefix} Removing {n_constants} constant BMKs (near-zero variance)")
            df = df[variable_mask].reset_index(drop=True).copy() 
            len_df = len(df)
            log.info(f"{log_prefix} {len_df} variable BMKs remaining")
            results["n_constants_removed"] = n_constants
        
        # Early exit if no variable BMKs
        if len_df == 0:
            log.warning(f"{log_prefix} No variable BMKs remaining after removing constants")
            return results
    elif stat_test == "beta-binomial":
        log.info(f"{log_prefix} Skipping proportion-variance prefilter for count-based inference")
    
    if is_active("filter2"):
        # PRE-FILTER #2: For very large sections, reduce to top 10k by variance
        # This prevents memory crashes in parallel workers for all_sequence_bmks / all_global_bmks
        MAX_BMKS_FOR_ANALYSIS = 10000
        if len_df > MAX_BMKS_FOR_ANALYSIS:
            log.info(f"{log_prefix} {len_df} BMKs exceeds {MAX_BMKS_FOR_ANALYSIS} limit")
            if stat_test == "beta-binomial":
                log.info(f"{log_prefix} Selecting top {MAX_BMKS_FOR_ANALYSIS} BMKs by total trials...")
                trial_cols = [f"{sample['col']}::trials" for sample in sample_info
                              if f"{sample['col']}::trials" in df.columns]
                coverage = df[trial_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
                top_indices = coverage.nlargest(MAX_BMKS_FOR_ANALYSIS).index
            else:
                log.info(f"{log_prefix} Selecting top {MAX_BMKS_FOR_ANALYSIS} BMKs by variance...")
                ndf_size_limit = numeric_df(df, sample_cols)
                variance_size_limit = ndf_size_limit[sample_cols].var(axis=1)
                top_indices = variance_size_limit.nlargest(MAX_BMKS_FOR_ANALYSIS).index
            df = df.loc[top_indices].reset_index(drop=True).copy()
            len_df = len(df)
            log.info(f"{log_prefix} Reduced to {len_df} BMKs")
            results["n_bmks_after_size_limit"] = len_df

        # Save filtered data
        df.to_csv(os.path.join(outdir, "data.csv"), index=True)
    

    # 1. Quality Control
    if is_active("qc"):
        try:
            log.info(f"{log_prefix} QC analysis...")
            results["qc"] = qc_analysis(df, sample_cols, os.path.join(outdir, "1_qc"))
        except Exception as e:
            log.warning(f"{log_prefix} QC failed: {e}")
    
    # 2. Batch Effect Detection
    if is_active("batch"):
        try:
            log.info(f"{log_prefix} Batch effect analysis...")
            results["batch"] = batch_effect_analysis(df, sample_cols, sample_info, os.path.join(outdir, "2_batch"))
        except Exception as e:
            log.warning(f"{log_prefix} Batch effect analysis failed: {e}")

    # 3. Descriptive Statistics
    if is_active("descriptive"):
        try:
            log.info(f"{log_prefix} Descriptive stats...")
            results["descriptive"] = descriptive_stats(df, sample_cols, sample_info, os.path.join(outdir, "3_descriptive"))
        except Exception as e:
            log.warning(f"{log_prefix} Descriptive stats failed: {e}")

    # 4. Multivariate Analysis (PCA, clustering)
    if is_active("multivariate"):
        try:
            log.info(f"{log_prefix} Multivariate analysis...")
            results["multivariate"] = multivariate_analysis(df, sample_cols, sample_info, os.path.join(outdir, "4_multivariate"))
        except Exception as e:
            log.warning(f"{log_prefix} Multivariate analysis failed: {e}")

    # 5. Differential Editing Analysis
    if is_active("differential"):
        try:
            log.info(f"{log_prefix} Differential analysis...")
            results["differential"] = differential_analysis(df, sample_cols, sample_info, os.path.join(outdir, "5_differential"), stat_test=stat_test)
        except Exception as e:
            if stat_test == "beta-binomial":
                raise
            log.warning(f"{log_prefix} Differential analysis failed: {e}")

    # 6. Correlation / Network Analysis
    if is_active("correlation"):
        try:
            log.info(f"{log_prefix} Correlation / Network analysis...")
            results["correlation"] = correlation_network(df, sample_cols, os.path.join(outdir, "6_correlation"))
        except Exception as e:
            log.warning(f"{log_prefix} Correlation analysis failed: {e}")

    # 7. Feature Selection / Biomarker Ranking
    if is_active("ranking"):
        try:
            log.info(f"{log_prefix} Feature ranking...")
            results["ranking"] = feature_ranking(df, sample_cols, sample_info, os.path.join(outdir, "7_ranking"), max_bmks=max_bmks, bmk_filter_cols=bmk_filter_cols)
        except Exception as e:
            log.warning(f"{log_prefix} Feature ranking failed: {e}")

    # 8. Classification / Predictive Modeling
    if is_active("classification"):
        try:
            log.info(f"{log_prefix} Classification...")
            results["classification"] = classification_analysis(df, sample_cols, sample_info, os.path.join(outdir, "8_classification"), max_bmks=max_bmks, bmk_filter_cols=bmk_filter_cols)
        except Exception as e:
            log.warning(f"{log_prefix} Classification failed: {e}")

    # 9. Stability / Robustness (replicate concordance)
    if is_active("stability"):
        try:
            log.info(f"{log_prefix} Stability analysis...")
            results["stability"] = stability_analysis(df, sample_cols, sample_info, os.path.join(outdir, "9_stability"))
        except Exception as e:
            log.warning(f"{log_prefix} Stability analysis failed: {e}")

    # 10. Heatmap
    if is_active("heatmap"):
        log.info(f"{log_prefix} Heatmap generation...")
        try:
            section_heatmap(df, sample_cols, sample_info, outdir, title=section_name)
            log.info(f"{log_prefix} Heatmap generation completed")
        except Exception as e:
            log.warning(f"{log_prefix} Heatmap failed: {e}")

    # END
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
