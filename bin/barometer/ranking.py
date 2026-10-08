"""barometer.ranking – Global biomarker ranking + significance diagrams.

global_ranking() aggregates the per-section differential tables into a single
global ranking (top-N plot, any/global/pairwise significant subsets).
generate_diagram() and split_significant_bmks() build the significance
heatmap and the unique/common per-comparison tree.
"""

import logging
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .utils import slurp_file, safe_mkdir, save_fig

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Global Biomarker Ranking (across all sections)
# ---------------------------------------------------------------------------

def global_ranking(all_results, outdir, stat_test="auto"):
    """Aggregate rankings across sections to produce a global ranking with cross-validation metrics."""
    safe_mkdir(outdir)

    # Collect differential results
    diff_files = []
    for vtype, vdata in all_results.items():
        for mtype, mdata in vdata.items():
            for section, sdata in mdata.items():
                # MEMORY FIX: sdata is now a slim dict with only "differential_table" key
                diff_path = sdata.get("differential_table")
                if diff_path and os.path.exists(diff_path):

                    ddf = slurp_file(diff_path, separator=",")
                    # MEMORY: no ddf.copy() — ddf is fresh from slurp_file, only 3 cols added
                    ddf["value_type"] = vtype
                    ddf["mtype"] = mtype
                    ddf["section"] = section
                    diff_files.append(ddf)

    if diff_files:
        all_diff = pd.concat(diff_files, ignore_index=True)
        build_global_ranking(all_diff, outdir, stat_test)

    return outdir


def build_global_ranking(all_diff, outdir, stat_test="auto"):
    """Build all global-ranking outputs from a concatenated differential DataFrame.

    ``all_diff`` must contain at least the columns: ``uid``, ``value_type``,
    ``mtype``, ``section``, ``primary_padj`` (plus the per-test ``*_pval`` /
    ``*_padj`` columns). This is the shared engine used both by
    :func:`global_ranking` (per mtype) and by the cross-mtype merge.
    """
    if all_diff is not None and len(all_diff) > 0:
        # Add n_sections_significant: cross-validation metric
        # Count how many sections each BMK appears in (regardless of significance)
        if "uid" in all_diff.columns:
            id_section_counts = all_diff.groupby("uid")["section"].nunique()
            all_diff["n_sections_total"] = all_diff["uid"].map(id_section_counts)
            
            # Count sections where BMK is significant (padj < 0.05 in ANY test)
            padj_cols = [c for c in all_diff.columns if c.endswith("_padj")]
            if padj_cols:
                # Create mask: True if ANY padj column is < 0.05
                sig_mask = all_diff[padj_cols].lt(0.05).any(axis=1)
                sig_df = all_diff[sig_mask].copy()
                if len(sig_df) > 0:
                    sig_counts = sig_df.groupby("uid")["section"].nunique()
                    all_diff["n_sections_significant"] = all_diff["uid"].map(sig_counts).fillna(0).astype(int)
                else:
                    all_diff["n_sections_significant"] = 0
            else:
                all_diff["n_sections_significant"] = 0

        # Rank by primary_padj
        if "primary_padj" in all_diff.columns:
            ranked = all_diff.dropna(subset=["primary_padj"]).sort_values("primary_padj")
            ranked.to_csv(os.path.join(outdir, "global_ranking_primary_padj.csv"), index=False)
            
            # Log statistics
            best_padj = ranked["primary_padj"].min() if len(ranked) > 0 else float('nan')
            log.info(f"  Saved global ranking: {len(ranked)} biomarkers (best padj: {best_padj:.2e})")

            # -----------------------------------------------------------------------
            # Top 50 plot (ranking complet) — reste au top-level
            # -----------------------------------------------------------------------
            top_n = min(50, len(ranked))
            top = ranked.head(top_n)
            fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.3)))
            neg_log_p = -np.log10(top["primary_padj"].clip(lower=1e-300))

            # Enhanced labels with n_sections_significant
            if "n_sections_significant" in top.columns:
                labels = [f"{row['uid']} [{row['section']}] (×{int(row['n_sections_significant'])} sec)" 
                          if row['n_sections_significant'] > 1 else f"{row['uid']} [{row['section']}]"
                          for _, row in top.iterrows()]
            else:
                labels = top["uid"].astype(str) + " [" + top["section"].astype(str) + "]"

            ax.barh(range(top_n), neg_log_p.values[::-1])
            ax.set_yticks(range(top_n))
            ax.set_yticklabels(labels[::-1], fontsize=6)
            ax.set_xlabel("-log10(adjusted p-value)")
            ax.set_title("Global BMK Ranking by Significance (top 50)\n(×N sec = significant in N sections)")

            # Add reference lines for p-value thresholds
            ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7, label='p=0.05')
            ax.axvline(-np.log10(0.1), color='orange', linestyle='--', linewidth=1, alpha=0.7, label='p=0.1')
            ax.legend(loc='lower right', fontsize=8)

            save_fig(fig, os.path.join(outdir, "global_ranking_plot.png"))
            plt.close(fig)

            # -----------------------------------------------------------------------
            # ANY comparison (union global ∪ pairwise)
            # → significant_bmks_any_comparison/
            # Consistent with the section-level "Found N significant biomarkers"
            # logic: significant if padj < 0.05 in ANY test (global or pairwise).
            # Both branches check the GLOBAL padj AND all PAIRWISE padjs:
            #   beta-binomial: beta_binomial_padj (global) + beta_binomial_padj_{pair} (pairwise)
            #   other tests:   primary_padj (global)       + primary_padj_{pair} (pairwise)
            # The only difference is the column prefix (beta_binomial_ vs primary_).
            # -----------------------------------------------------------------------
            any_dir = os.path.join(outdir, "significant_bmks_any_comparison")
            safe_mkdir(any_dir)

            if stat_test == "beta-binomial":
                sig_padj_cols = [c for c in ranked.columns
                                 if c == "beta_binomial_padj" or c.startswith("beta_binomial_padj_")]
                sig_mask = ranked[sig_padj_cols].lt(0.05).any(axis=1)
                sig_ranked = ranked[sig_mask].copy()
                sig_desc = "padj < 0.05 in any beta-binomial test (global or pairwise)"
            else:
                sig_padj_cols = [c for c in ranked.columns
                                 if c == "primary_padj" or c.startswith("primary_padj_")]
                sig_mask = ranked[sig_padj_cols].lt(0.05).any(axis=1)
                sig_ranked = ranked[sig_mask].copy()
                sig_desc = "padj < 0.05 in any test (global or pairwise)"

            if len(sig_ranked) > 0:
                sig_csv_path = os.path.join(any_dir, "significant_bmks.csv")
                sig_ranked.to_csv(sig_csv_path, index=False)
                log.info(f"  Saved any-comparison subset: {len(sig_ranked)} biomarkers ({sig_desc})")
                # -----------------------------------------------------------------------
                # Create figure from significant_bmks.csv
                # -----------------------------------------------------------------------
                try:
                    generate_diagram(sig_csv_path, any_dir)
                except Exception as e:
                    log.warning(f"  generate_diagram failed: {e}")
                # -----------------------------------------------------------------------
                # Create structure folder with unique/common splits for each comparison
                # -----------------------------------------------------------------------
                try:
                    split_significant_bmks(sig_csv_path, os.path.join(outdir, "significant_bmks_pairwise_comparison"))
                except Exception as e:
                    log.warning(f"  split_significant_bmks failed: {e}")
            else:
                log.info(f"  No significant biomarkers found ({sig_desc})")

            # -----------------------------------------------------------------------
            # GLOBAL test only (omnibus)
            # → significant_bmks_global_comparison/
            # primary_padj is ALWAYS the global test padj (Kruskal/ANOVA/Welch in
            # normal modes, beta_binomial_padj in beta-binomial mode), so this is
            # mode-agnostic: significant if the global test padj < 0.05.
            # -----------------------------------------------------------------------
            global_dir = os.path.join(outdir, "significant_bmks_global_comparison")
            safe_mkdir(global_dir)

            global_sig = ranked[ranked["primary_padj"] < 0.05].copy()
            if len(global_sig) > 0:
                global_sig_csv = os.path.join(global_dir, "significant_bmks.csv")
                global_sig.to_csv(global_sig_csv, index=False)
                log.info(f"  Saved global-significant subset: {len(global_sig)} biomarkers (primary_padj < 0.05)")

                # Heatmap de significativité (complément du bar chart ci-dessous)
                try:
                    generate_diagram(global_sig_csv, global_dir)
                except Exception as e:
                    log.warning(f"  generate_diagram (global) failed: {e}")

                # Top 50 plot for global-significant BMKs
                top_n_g = min(50, len(global_sig))
                top_g = global_sig.head(top_n_g)
                fig_g, ax_g = plt.subplots(figsize=(8, max(4, top_n_g * 0.3)))
                neg_log_p_g = -np.log10(top_g["primary_padj"].clip(lower=1e-300))
                if "n_sections_significant" in top_g.columns:
                    labels_g = [f"{row['uid']} [{row['section']}] (×{int(row['n_sections_significant'])} sec)" 
                                if row['n_sections_significant'] > 1 else f"{row['uid']} [{row['section']}]"
                                for _, row in top_g.iterrows()]
                else:
                    labels_g = top_g["uid"].astype(str) + " [" + top_g["section"].astype(str) + "]"
                ax_g.barh(range(top_n_g), neg_log_p_g.values[::-1])
                ax_g.set_yticks(range(top_n_g))
                ax_g.set_yticklabels(labels_g[::-1], fontsize=6)
                ax_g.set_xlabel("-log10(adjusted p-value)")
                ax_g.set_title("Global-Significant BMKs (top 50)\n(primary_padj < 0.05)")
                ax_g.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7, label='p=0.05')
                ax_g.axvline(-np.log10(0.1), color='orange', linestyle='--', linewidth=1, alpha=0.7, label='p=0.1')
                ax_g.legend(loc='lower right', fontsize=8)
                save_fig(fig_g, os.path.join(global_dir, "global_ranking_plot.png"))
                plt.close(fig_g)
            else:
                log.info(f"  No global-significant biomarkers (primary_padj < 0.05)")

# ---------------------------------------------------------------------------
# Venn Diagram Generation
# ---------------------------------------------------------------------------

def generate_diagram(sig_csv_path, outdir):
    """Génère un diagramme de Venn pour les sections présentes dans global_ranking_significant.csv"""
    safe_mkdir(outdir)

    # Charger le fichier
    df = pd.read_csv(sig_csv_path)
    log.debug(f"Loaded {len(df)} rows from {sig_csv_path}")
    log.debug(f"Columns: {df.columns.tolist()}")

    # 🔎 1. sélectionner les colonnes de padj (comparaisons)
    padj_cols = [c for c in df.columns if c.startswith("primary_padj_")]
    log.debug(f"Found {len(padj_cols)} padj columns: {padj_cols}")

    if len(padj_cols) == 0:
        log.warning("No primary_padj columns found")
        return

    # 🧼 2. nettoyer les noms de colonnes
    rename_map = {
        c: c.replace("primary_padj_", "")
        for c in padj_cols
    }
    df = df.rename(columns=rename_map)
    padj_cols = list(rename_map.values())
    log.debug(f"Renamed columns: {padj_cols}")

    # 🎯 3. Créer une matrice avec 3 niveaux de significativité
    sig_matrix = pd.DataFrame(0, index=df.index, columns=padj_cols)
    
    for col in padj_cols:
        pvals = df[col]
        log.debug(f"{col} - min={pvals.min():.4f}, max={pvals.max():.4f}, NaN={pvals.isna().sum()}")
        sig_matrix.loc[pvals < 0.1, col] = 1   # marginalement significatif
        sig_matrix.loc[pvals < 0.05, col] = 2  # significatif
    
    log.debug(f"sig_matrix value counts:\n{sig_matrix.stack().value_counts()}")
    
    # Utiliser uid comme index si disponible
    if "uid" in df.columns:
        sig_matrix.index = df["uid"]
        log.debug(f"Using uid as index")
    
    # 📊 4. plot avec palette à 3 couleurs
    n_rows = sig_matrix.shape[0]
    fig_height = max(6, n_rows * 0.3)
    log.debug(f"Creating figure {12}x{fig_height} for {n_rows} rows")
    
    fig, ax = plt.subplots(figsize=(12, fig_height))
    
    # Import nécessaire
    from matplotlib.colors import ListedColormap
    
    # Palette personnalisée : blanc (non-sig), orange (p<0.1), rouge (p<0.05)
    colors = ["#f0f0f0", "#FFA500", "#D32F2F"]
    cmap = ListedColormap(colors)
    
    # Heatmap
    sns.heatmap(
        sig_matrix, 
        cmap=cmap, 
        ax=ax,
        vmin=0,
        vmax=2,
        linewidths=0.5,
        linecolor='lightgray',
         cbar=False 
      #  cbar_kws={
      #      'label': 'Significance',
      #      'ticks': [0.33, 1, 1.67],  # ⭐ Position des ticks au centre de chaque couleur
      #      'boundaries': [0, 0.66, 1.33, 2]  # ⭐ Limites entre couleurs
      #  }
    )
    
    # colorbar label
    #cbar = ax.collections[0].colorbar
    #cbar.set_ticklabels(['NS\n(p≥0.1)', 'Marginal\n(p<0.1)', 'Significant\n(p<0.05)'])

    # Labels
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0, fontsize=8)
    plt.xlabel("Comparisons", fontsize=10)
    plt.ylabel("Biomarkers (uid)", fontsize=10)
    plt.title("Significance of biomarkers across conditions", fontsize=12, pad=20)
    
    # Légende manuelle
    from matplotlib.patches import Rectangle
    legend_elements = [
        Rectangle((0, 0), 1, 1, fc=colors[0], label='Not significant (p ≥ 0.1)'),
        Rectangle((0, 0), 1, 1, fc=colors[1], label='Marginal (0.05 ≤ p < 0.1)'),
        Rectangle((0, 0), 1, 1, fc=colors[2], label='Significant (p < 0.05)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.15, 1))

    # 💾 Save
    output_path = os.path.join(outdir, "ranking_plot.png")
    log.debug(f"Saving to {output_path}")
    plt.savefig(
        output_path, 
        bbox_inches="tight", 
        facecolor="white",
        dpi=150
    )
    plt.close()
    
    log.info(f"✓ Saved significance heatmap with {n_rows} BMKs × {len(padj_cols)} comparisons")
    log.info(f"  Saved significance heatmap: {output_path}")

# ---------------------------------------------------------------------------
# Split csv significant BMKs into unique/common for each comparison
# ---------------------------------------------------------------------------
def split_significant_bmks(sig_csv_path, outdir):
    """Crée une arborescence de dossiers/fichiers pour chaque comparaison :
    - all_significant.csv : tous les BMKs significatifs pour la comparaison
    - unique/unique.csv : BMKs uniquement significatifs pour cette comparaison
    - common/<autre_comparaison>.csv : BMKs partagés avec une autre comparaison
    """
    safe_mkdir(outdir)

    df = pd.read_csv(sig_csv_path)
    df_orig = df.copy()  # garde les colonnes primary_padj_* pour le plot
    padj_cols = [c for c in df.columns if c.startswith("primary_padj_")]
    rename_map = {c: c.replace("primary_padj_", "") for c in padj_cols}
    df = df.rename(columns=rename_map)
    padj_cols = list(rename_map.values())
    alpha = 0.05
    binary_df = df[padj_cols] < alpha
    binary_df = binary_df.set_index(df["uid"])  # index = BMK ID

    for comp in padj_cols:
        comp_dir = os.path.join(outdir, comp)
        # 1. Tous les BMKs significatifs pour cette comparaison
        sig_bmks = binary_df.index[binary_df[comp]].tolist()
        all_sig_df = df[df[comp] < alpha]
        if not all_sig_df.empty:
            os.makedirs(comp_dir, exist_ok=True)
            all_sig_df.to_csv(os.path.join(comp_dir, "all_significant.csv"), index=False)
            # ranking_plot.png : heatmap de significativité des BMKs de cette comparaison
            orig_col = rename_map[comp]  # ex: "primary_padj_long.non.cleared_vs_negative"
            comp_sig_orig = df_orig[df_orig[orig_col] < alpha]
            if not comp_sig_orig.empty:
                comp_sig_csv = os.path.join(comp_dir, "_ranking_input.csv")
                comp_sig_orig.to_csv(comp_sig_csv, index=False)
                try:
                    generate_diagram(comp_sig_csv, comp_dir)
                    os.remove(comp_sig_csv)  # on ne garde que le plot
                except Exception as e:
                    log.warning(f"  generate_diagram ({comp}) failed: {e}")
        # 2. BMKs uniques à cette comparaison
        is_unique = (binary_df[comp]) & (binary_df.drop(columns=[comp]).sum(axis=1) == 0)
        unique_df = df[df["uid"].isin(binary_df.index[is_unique])]
        if not unique_df.empty:
            unique_dir = os.path.join(comp_dir, "unique")
            os.makedirs(unique_dir, exist_ok=True)
            unique_df.to_csv(os.path.join(unique_dir, "unique.csv"), index=False)
        # 3. BMKs communs avec chaque autre comparaison
        for other in padj_cols:
            if other == comp:
                continue
            is_common = (binary_df[comp]) & (binary_df[other])
            common_df = df[df["uid"].isin(binary_df.index[is_common])]
            if not common_df.empty:
                common_dir = os.path.join(comp_dir, "common")
                os.makedirs(common_dir, exist_ok=True)
                common_df.to_csv(os.path.join(common_dir, f"{other}.csv"), index=False)
