"""barometer.__main__ – Entry point for `python -m barometer`.

Contains the CLI (main) and the picklable analyze_section_wrapper used by
the ProcessPoolExecutor.
"""

import argparse
import gc
import glob
import json
import logging
import multiprocessing
import os
import shutil
import sys
import time
import warnings

# Limit implicit BLAS/LAPACK multi-threading so each of the -j worker
# processes uses at most 1 core (total CPU usage = n_jobs, not n_jobs^2).
# Must be set before numpy/scipy/sklearn are imported.
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')

from concurrent.futures import ProcessPoolExecutor, as_completed
# The container runs Python 3.12, where concurrent.futures.TimeoutError IS the
# builtin TimeoutError (unified in 3.11). A plain `except TimeoutError:` is
# sufficient; no separate futures import is needed.

try:
    import psutil
except ImportError:
    psutil = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .analysis import analyze_section, build_feature_tree
from .ranking import global_ranking
from .utils import (
    slurp_file,
    safe_mkdir,
    harmonize_columns,
    create_uid_column,
    parse_sample_columns,
    get_value_types,
    cols_for_vtype,
    sample_info_for_vtype,
    prepare_df_for_task,
    filter_significant,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def analyze_section_wrapper(args_tuple):
    """Wrapper for analyze_section to be used with ProcessPoolExecutor.
    
    Args:
        args_tuple: (df_source, sample_cols, sample_info, outdir, section_name, result_key, stat_test, bmk_filter_cols, max_bmks, task_timeout)
            where df_source can be:
            - dict: legacy format, converted back to DataFrame
            - tuple: ('pickled', bytes_data) for memory-efficient transfer
            task_timeout: per-task SIGALRM budget in seconds (doubled on retry)
        
    Returns:
        (result_key, results_dict)
    """
    import traceback
    import signal
    import gc
    import pickle
    import random
    import time as time_module
    df_source, sample_cols, sample_info, outdir, section_name, result_key, stat_test, bmk_filter_cols, max_bmks, task_timeout = args_tuple
    
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
    
    # Setup timeout alarm (per-task budget, doubled by the supervisor on retry)
    def timeout_handler(signum, frame):
        raise TimeoutError(f"Task exceeded {task_timeout} seconds: {section_name}")
    
    try:
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(task_timeout)
        
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
        # First pass now runs ALL steps per section (QC/batch/descriptive on
        # all BMKs, then differential, then selection-dependent steps on the
        # differential BMKs — see analyze_section).
        results = analyze_section(df, sample_cols, sample_info, outdir, section_name, stat_test=stat_test, bmk_filter_cols=bmk_filter_cols, max_bmks=max_bmks, enabled_tests=None)
        
        # Cancel alarm
        signal.alarm(0)
        
        # MEMORY FIX #2: Return slim results dict (only file paths, not full data)
        # Full results dict is ~500KB per task × 1491 tasks = 745MB accumulated in all_results
        # global_ranking() only needs the differential.table path (50 bytes)
        slim_results = {
            "differential_table": os.path.join(outdir, "5_differential", "differential_results.csv")
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
  
  BETA-BINOMIAL MODE (--stat-test beta-binomial)
    Instead of testing the rounded proportions, this mode models the raw,
    unrounded DRIP counts per sample:
        {sample}::successes  = reads showing the edited base
        {sample}::trials     = total reads covering the site
    (These count columns must be present, i.e. DRIP must have been run with
    count columns enabled. Requires Rscript + glmmTMB → use the Barometer
    container.)
    
    For each biomarker a single beta-binomial GLMM is fitted (R/glmmTMB):
        cbind(successes, trials - successes) ~ condition
    The beta-binomial dispersion parameter absorbs the extra-binomial
    overdispersion typical of sequencing data, giving more reliable p-values
    than a plain binomial or a test on rounded proportions.
    
    Two kinds of Wald tests are derived from that one model:
      • GLOBAL   (1 per BMK): "does condition have ANY effect?"
          → beta_binomial_pval / beta_binomial_padj   (like an ANOVA F-test)
      • PAIRWISE (C(n,2) per BMK): one contrast per group pair
          → beta_binomial_pval_{g1}_vs_{g2} / beta_binomial_padj_{g1}_vs_{g2}
    
    A BMK can be significant in one pairwise contrast while the global test is
    not (effect concentrated on a single pair). Significance = padj < 0.05 in
    ANY test (global or pairwise). In this mode primary_* columns are populated
    from the beta-binomial results and the near-zero-variance prefilter is
    skipped.
  
  OPTION 1: --stat-test (controls which test is emphasized in plots)
    This creates primary_pval and primary_padj columns that COPY the selected test:
    
    nonparametric : primary_* = kruskal_* (DEFAULT, robust)
    parametric    : primary_* = anova_*
    welch         : primary_* = welch_*
    auto          : primary_* = auto-selected test (varies per BMK)
    kruskal       : alias for nonparametric
    beta-binomial : primary_* = beta_binomial_* (see BETA-BINOMIAL MODE above)
  
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
  ./barometer_analyze.py -a data.tsv -o results/ \\
    --stat-test welch --bmk-filter welch_padj kruskal_padj --max-bmks 1000
  
  # Conservative: only most significant BMKs (200 max):
  ./barometer_analyze.py -a data.tsv -o results/ --max-bmks 200
  
  # Aggressive: maximize BMK collection with high limit:
  ./barometer_analyze.py -a data.tsv -o results/ --max-bmks 1500
  
  # Beta-binomial on raw DRIP counts (needs Rscript + glmmTMB, Barometer container):
  ./barometer_analyze.py -a data.tsv -o results/ --stat-test beta-binomial
        """)
    parser.add_argument("-a", "--aggregates", default=None, help="Aggregates TSV file (optional)")
    parser.add_argument("-f", "--features", default=None, help="Features TSV file (optional)")
    parser.add_argument("-o", "--outdir", default="barometer_results", help="Output directory")
    parser.add_argument("--sites", default=None, help="Per-site DRIP TSV (drip.py on the per-site pluviometer output); analysed as a single section, in place of --features")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Number of parallel jobs (default: 1, use -1 for all CPUs)")
    parser.add_argument("-v", "--value-types", nargs="+", default=None, help="Value types to analyze (e.g., espf espr). If not specified, all value types are analyzed.")
    parser.add_argument("--agg-levels", nargs="+", default=None, help="Aggregate levels to analyze: global, sequence, feature. If not specified, all levels are analyzed.")
    parser.add_argument("--feature-types", nargs="+", default=None, help="Feature types to analyze (e.g., gene exon RNA). If not specified, all feature types are analyzed.")
    parser.add_argument("--stat-test", default="nonparametric", 
                        choices=["auto", "parametric", "nonparametric", "welch", "kruskal", "beta-binomial"],
                        help="Statistical test selection: beta-binomial (requires DRIP counts), auto, parametric, nonparametric, welch or kruskal")
    parser.add_argument("--bmk-filter", nargs="+", default=None,
                        help="Priority list for significant BMK selection (default: primary_padj only for beta-binomial; legacy cascade otherwise).")
    parser.add_argument("--max-bmks", type=int, default=500,
                        help="Maximum number of BMKs to use for RandomForest and classification (default: 500). Prevents memory crashes with large datasets.")
    parser.add_argument("--merge", action="store_true",
                        help="Merge mode: combine the per-mtype global rankings into a truly-global ranking and run a second pass by redistribution per value type.")
    parser.add_argument("--results-dir", nargs="+", default=None,
                        help="In --merge mode: one or more barometer_results dirs, each containing "
                             "espf/ and/or espr/ subdirs with per-mtype (aggregate/feature/sites) global_ranking CSVs.")
    parser.add_argument("--raw-inputs", default=None,
                        help="In --merge mode: JSON mapping of vtype -> raw TSVs, e.g. "
                             '\'{"espf": {"aggregates": "a.tsv", "features": "f.tsv", "sites": "s.tsv"}, '
                             '"espr": {"aggregates": "a.tsv", "features": "f.tsv", "sites": "s.tsv"}}\'. '
                             "Used to rebuild full rows for the second pass.")
    parser.add_argument("--resume", default=None,
                        help="Resume mode: path to a failed_analyses.tsv from a previous run. "
                             "Only the listed sections (vtype + mtype + section_key) are re-run. "
                             "The new failed_analyses.tsv then contains only the sections that "
                             "still failed, so --resume can be iterated until it is empty.")
    args = parser.parse_args()

    # Resume mode: load the set of (vtype, mtype, section_key) to re-run.
    resume_keys = None
    if args.resume:
        if not os.path.isfile(args.resume):
            log.error(f"--resume file not found: {args.resume}")
            sys.exit(1)
        df_failed = pd.read_csv(args.resume, sep="\t", low_memory=False)
        required_cols = {"vtype", "mtype", "section_key"}
        missing = required_cols - set(df_failed.columns)
        if missing:
            log.error(f"--resume file {args.resume} is missing columns: {sorted(missing)}")
            sys.exit(1)
        resume_keys = set(zip(df_failed["vtype"].astype(str),
                              df_failed["mtype"].astype(str),
                              df_failed["section_key"].astype(str)))
        log.info(f"RESUME MODE: {len(resume_keys)} sections to re-run from {args.resume}")

    # ------------------------------------------------------------------
    # MERGE MODE: combine per-mtype global rankings into a truly-global one.
    # ------------------------------------------------------------------
    if args.merge:
        if not args.results_dir:
            parser.error("--merge requires --results-dir")
        if args.bmk_filter is None:
            args.bmk_filter = (
                ["primary_padj"]
                if args.stat_test == "beta-binomial"
                else ["primary_padj", "kruskal_padj", "welch_padj", "anova_padj"]
            )
        raw_inputs = None
        if args.raw_inputs:
            raw_inputs = json.loads(args.raw_inputs)
        from .merge import merge_global_rankings
        merge_global_rankings(
            results_dirs=args.results_dir,
            outdir=args.outdir,
            stat_test=args.stat_test,
            raw_inputs=raw_inputs,
            max_bmks=args.max_bmks,
            bmk_filter=args.bmk_filter,
        )
        return

    if args.stat_test == "beta-binomial" and shutil.which("Rscript") is None:
        parser.error("Beta-binomial analysis requires Rscript and glmmTMB; use the Barometer container.")

    if args.bmk_filter is None:
        args.bmk_filter = (
            ["primary_padj"]
            if args.stat_test == "beta-binomial"
            else ["primary_padj", "kruskal_padj", "welch_padj", "anova_padj"]
        )
    
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
    if not args.aggregates and not args.features and not args.sites:
        parser.error("At least one of --aggregates, --features or --sites must be provided")

    outdir = args.outdir
    safe_mkdir(outdir)

    # Load data
    log.info("Loading data...")
    agg_df = None
    feat_df = None
    sites_df = None
    # mtypes whose file was provided but had no data rows: no analysis is
    # run, but the mtype dir + manifest are still written so downstream
    # consumers (Nextflow) find the expected output.
    empty_mtypes = set()

    if args.aggregates:
        if not os.path.exists(args.aggregates):
            log.error(f"Aggregates file not found: {args.aggregates}")
            sys.exit(1)
        log.info(f"Loading aggregates from {args.aggregates}...")
        agg_df = slurp_file(args.aggregates, separator="\t")
        agg_df.columns = agg_df.columns.str.strip()  # Remove leading/trailing whitespace from column names
        agg_df = harmonize_columns(agg_df)  # Ensure column header standardisation
        agg_df = create_uid_column(agg_df)  # Create UID column for unique identification of rows

        # Strip whitespace from string columns
        for col in agg_df.select_dtypes(include=['object', 'string']).columns:
            agg_df[col] = agg_df[col].str.strip() if agg_df[col].dtype in ['object', 'string'] else agg_df[col]
        # MEMORY: convert low-cardinality columns to category (saves ~50-90 bytes/row/col)
        for col in ("Mtype", "Ptype", "Type", "Ctype", "Mode", "Strand", "SeqID"):
            if col in agg_df.columns:
                agg_df[col] = agg_df[col].astype("category")
        log.info(f"  Aggregates: {len(agg_df)} rows")
        if len(agg_df) == 0:
            log.warning("  Aggregates file has no data rows — no analysis will be run "
                        "(mtype dir + manifest are still written).")
            empty_mtypes.add("aggregate")
    else:
        log.info("Skipping aggregates (no file provided)")
    
    if args.features:
        if not os.path.exists(args.features):
            log.error(f"Features file not found: {args.features}")
            sys.exit(1)
        log.info(f"Loading features from {args.features}...")
        feat_df = slurp_file(args.features, separator="\t") 
        feat_df.columns = feat_df.columns.str.strip()  # Remove leading/trailing whitespace from column names
        feat_df = harmonize_columns(feat_df)  # Ensure column header standardisation
        feat_df = create_uid_column(feat_df)  # Create UID column for unique identification of rows
        # Strip whitespace from string columns
        for col in feat_df.select_dtypes(include=['object', 'string']).columns:
            feat_df[col] = feat_df[col].str.strip() if feat_df[col].dtype in ['object', 'string'] else feat_df[col]
        # MEMORY: convert low-cardinality columns to category (saves ~50-90 bytes/row/col)
        for col in ("Mtype", "Ptype", "Type", "Ctype", "Mode", "Strand", "SeqID"):
            if col in feat_df.columns:
                feat_df[col] = feat_df[col].astype("category")
        log.info(f"  Features: {len(feat_df)} rows")
        if len(feat_df) == 0:
            log.warning("  Features file has no data rows — no analysis will be run "
                        "(mtype dir + manifest are still written).")
            empty_mtypes.add("feature")
    else:
        log.info("Skipping features (no file provided)")

    if args.sites:
        if not os.path.exists(args.sites):
            log.error(f"Sites file not found: {args.sites}")
            sys.exit(1)
        log.info(f"Loading sites from {args.sites}...")
        sites_df = slurp_file(args.sites, separator="\t")
        sites_df.columns = sites_df.columns.str.strip()  # Remove leading/trailing whitespace from column names
        sites_df = harmonize_columns(sites_df)  # Ensure column header standardisation
        sites_df = create_uid_column(sites_df)  # Create UID column for unique identification of rows
        # Strip whitespace from string columns
        for col in sites_df.select_dtypes(include=['object', 'string']).columns:
            sites_df[col] = sites_df[col].str.strip() if sites_df[col].dtype in ['object', 'string'] else sites_df[col]
        # MEMORY: convert low-cardinality columns to category (saves ~50-90 bytes/row/col)
        for col in ("Mtype", "Ptype", "Type", "Ctype", "Mode", "Strand", "SeqID"):
            if col in sites_df.columns:
                sites_df[col] = sites_df[col].astype("category")
        log.info(f"  Sites: {len(sites_df)} rows")
        if len(sites_df) == 0:
            log.warning("  Sites file has no data rows — no analysis will be run "
                        "(mtype dir + manifest are still written).")
            empty_mtypes.add("sites")
    else:
        log.info("Skipping sites (no file provided)")

    # Parse sample information from whichever file is available
    sample_df = agg_df if agg_df is not None else (feat_df if feat_df is not None else sites_df)
    sample_info = parse_sample_columns(sample_df.columns)
    if args.stat_test == "beta-binomial":
        missing_counts = []
        for data_frame in (agg_df, feat_df, sites_df):
            if data_frame is None:
                continue
            frame_samples = parse_sample_columns(data_frame.columns)
            missing_counts.extend(
                f"{sample['col']}::{suffix}"
                for sample in frame_samples
                for suffix in ("successes", "trials")
                if f"{sample['col']}::{suffix}" not in data_frame.columns
            )
        if missing_counts:
            parser.error(
                "Beta-binomial analysis requires count columns in DRIP TSVs; "
                f"missing examples: {', '.join(missing_counts[:4])}. Rerun DRIP."
            )
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
    # Accumulates every failed/cancelled section so the user can re-run them
    # later. Written to <outdir>/failed_analyses.tsv at the end of the run.
    failed_analyses = []

    for vtype in value_types:
        log.info(f"\n{'='*60}")
        log.info(f"VALUE TYPE: {vtype}")
        log.info(f"{'='*60}")

        # Results land in <outdir>/<vtype>/<mtype>/ (mtype in aggregate/feature/sites).
        # In site mode the single section is stored under the "sites" mtype, so the
        # per-site results live alongside the aggregates/features of the same vtype.
        vtype_dir = os.path.join(outdir, vtype)
        safe_mkdir(vtype_dir)
        vcols = cols_for_vtype(sample_info, vtype)
        v_sample_info = sample_info_for_vtype(sample_info, vtype)
        all_results[vtype] = {"aggregate": {}, "feature": {}, "sites": {}}
        
        # Collect all analysis tasks for this value_type
        tasks = []

        # ===============================================================
        # AGGREGATES
        # ===============================================================
        if agg_df is not None and len(agg_df) > 0:
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
                    glob_agg, vcols, v_sample_info,
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
                    section, vcols, v_sample_info,
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
                        section, vcols, v_sample_info,
                        os.path.join(agg_dir, "global", f"-{mode}"),
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
                            section, vcols, v_sample_info,
                            os.path.join(agg_dir, "global", f"{ptype}-{mode}"),
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
                        section, vcols, v_sample_info,
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
                            section, vcols, v_sample_info,
                            os.path.join(agg_dir, f"sequence/{chrom}", f"{mode}"),
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
                                section, vcols, v_sample_info,
                                os.path.join(agg_dir, f"sequence/{chrom}", f"{ptype}-{mode}"),
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
                    section, vcols, v_sample_info,
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
                        section, vcols, v_sample_info,
                        os.path.join(agg_dir, "sequence", "all_sequence_bmks", f"{mode}"),
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
                            section, vcols, v_sample_info,
                            os.path.join(agg_dir, "sequence", "all_sequence_bmks", f"{ptype}-{mode}"),
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
                                section, vcols, v_sample_info,
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
                                    section, vcols, v_sample_info,
                                    os.path.join(agg_dir, "feature", "by_sequence", str(chrom), safe_name),
                                    f"Feature Agg - Chr {chrom} - Ptype={ptype}, Ctype={ctype}, Mode={mode}",
                                    ("aggregate", f"featagg_chr{chrom}_{safe_name}"),
                                    args.stat_test,
                                    args.bmk_filter,
                                    args.max_bmks
                                ))
        else:
            log.info(f"\n--- Skipping AGGREGATES for {vtype} (no aggregates data) ---")

        # ===============================================================
        # FEATURES
        # ===============================================================
        if feat_df is not None and len(feat_df) > 0:
            log.info(f"\n--- FEATURES for {vtype} ---")
            feat_data = feat_df  # Mtype is always "feature" in features file (no copy: read-only usage)
            feat_dir = os.path.join(vtype_dir, "feature")

            # Build hierarchy
            tree, top_ids = build_feature_tree(feat_data)
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
                    type_all_df, vcols, v_sample_info,
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
                        sub_df, vcols, v_sample_info,
                        os.path.join(type_dir, safe_fid),
                        f"Feature: {fid}",
                        ("feature", f"feature_{safe_fid}"),
                        args.stat_test,
                        args.bmk_filter,
                        args.max_bmks
                    ))
        else:
            log.info(f"\n--- Skipping FEATURES for {vtype} (no features data) ---")

        # ===============================================================
        # SITES
        # ===============================================================
        if sites_df is not None and len(sites_df) > 0:
            log.info(f"\n--- SITES for {vtype} ---")
            sites_data = sites_df
            tasks.append((
                sites_data, vcols, v_sample_info,
                os.path.join(vtype_dir, "sites", "all_sites"),
                "Sites - all",
                ("sites", "sites_all"),
                args.stat_test,
                args.bmk_filter,
                args.max_bmks
            ))
        else:
            log.info(f"\n--- Skipping SITES for {vtype} (no sites data) ---")

        # Resume mode: keep only the sections listed in failed_analyses.tsv.
        if resume_keys is not None:
            n_before = len(tasks)
            tasks = [t for t in tasks if (vtype, str(t[5][0]), str(t[5][1])) in resume_keys]
            log.info(f"  RESUME: {len(tasks)}/{n_before} sections selected for {vtype}")
            if not tasks:
                log.info(f"  RESUME: nothing to re-run for {vtype}, skipping")
                continue

        # MEMORY: tasks hold raw section DataFrames (views of the already-held
        # frames, no extra data). Pickling is deferred to submission time so
        # only max_pending_tasks pickled blobs exist in memory at once.
        if psutil is not None:
            try:
                log.info(f"  Parent RSS after task building: {psutil.Process().memory_info().rss / (1024 * 1024):.0f}MB ({len(tasks)} tasks)")
            except Exception:
                pass

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
            
            def record_failure(section_name, key, section_outdir, reason, detail=""):
                """Append a failed/cancelled section to the resume log."""
                mtype, section_key = key
                failed_analyses.append({
                    "vtype": vtype,
                    "mtype": mtype,
                    "section_key": section_key,
                    "section_name": section_name,
                    "section_outdir": section_outdir,
                    "reason": reason,
                    "detail": detail,
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                })

            # --- Retry policy -------------------------------------------------
            # A section that fails for a TIME-related reason (timeout, slow,
            # deadlock, orphaned) is retried up to MAX_RETRIES times, doubling
            # its per-task timeout each time: 300s -> 600s -> 1200s -> 2400s.
            # Code errors (error, crash) are NOT retried — more time won't fix
            # a bug. Only definitive failures (after retries exhausted) are
            # written to failed_analyses.tsv.
            BASE_TASK_TIMEOUT = 300
            MAX_RETRIES = 3
            TIME_RELATED_REASONS = {"timeout", "slow", "deadlock", "orphaned"}

            with ProcessPoolExecutor(max_workers=n_jobs, max_tasks_per_child=20, mp_context=mp_context) as executor:
                future_to_key = {}
                task_iter = iter(tasks)
                retry_queue = []          # 9-tuples re-submitted with a doubled timeout
                retry_count = {}          # key -> number of retries already used
                task_timeout_by_key = {}  # key -> current per-task timeout (doubled on retry)
                completed = 0
                failed = 0
                submitted_count = 0
                last_progress_time = time.time()
                last_worker_check = time.time()
                worker_check_interval = 15  # Check worker health every 15 seconds
                # Track submission time per future for per-future timeout
                submit_time = {}

                def current_future_timeout():
                    # Per-future timeout = the LARGEST per-task timeout currently
                    # in flight + grace for spawn/pickle overhead. A single slow
                    # task must NOT cancel the whole batch — only the future that
                    # actually exceeded its own budget is cancelled.
                    in_flight = [task_timeout_by_key.get(k, BASE_TASK_TIMEOUT) for k, _, _ in future_to_key.values()]
                    return (max(in_flight) if in_flight else BASE_TASK_TIMEOUT) + 30

                def current_stall_timeout():
                    # Global stall safety net: only fires if NO future completed
                    # for longer than the per-future timeout (true deadlock, e.g.
                    # all workers stuck). Must be > future_timeout.
                    return current_future_timeout() + 60

                def submit_task(task_tuple):
                    """Submit one 9-tuple task, applying its (possibly doubled) timeout."""
                    nonlocal submitted_count
                    df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks = task_tuple
                    tto = task_timeout_by_key.get(key, BASE_TASK_TIMEOUT)
                    # Pickle at submission time (not up-front) to keep
                    # parent memory bounded to max_pending_tasks blobs.
                    future = executor.submit(analyze_section_wrapper, (prepare_df_for_task(df_source), cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks, tto))
                    future_to_key[future] = (key, name, outdir)
                    submit_time[future] = time.time()
                    pending_futures.add(future)
                    submitted_count += 1
                    # Explicitly free the df_source reference to help GC
                    del df_source
                    return future

                def next_task():
                    """Pop the next task to submit: retries first, then the main queue."""
                    if retry_queue:
                        return retry_queue.pop(0)
                    return next(task_iter)

                def handle_time_failure(key, section_name, section_outdir, reason, detail):
                    """Retry a time-related failure (doubling the timeout) or record it.

                    Returns True if the task was re-queued for retry, False if it
                    was recorded as a definitive failure.
                    """
                    nonlocal failed
                    n = retry_count.get(key, 0)
                    if n < MAX_RETRIES:
                        retry_count[key] = n + 1
                        task_timeout_by_key[key] = BASE_TASK_TIMEOUT * (2 ** (n + 1))
                        # Rebuild the 9-tuple from the original task list by key.
                        retry_queue.append(task_by_key[key])
                        log.warning(f"  ↻ RETRY ({n+1}/{MAX_RETRIES}) with {task_timeout_by_key[key]}s timeout: {section_name} ({reason})")
                        return True
                    failed += 1
                    record_failure(section_name, key, section_outdir, reason, detail)
                    return False

                # Index tasks by key so a retry can re-fetch its original 9-tuple
                # (the df view is still alive in the parent; pickling is idempotent).
                task_by_key = {t[5]: t for t in tasks}

                # Process futures with stall detection
                pending_futures = set()

                def refill():
                    """Submit tasks until max_pending_tasks are in flight (or queues empty)."""
                    while len(pending_futures) < max_pending_tasks:
                        try:
                            submit_task(next_task())
                        except StopIteration:
                            break

                # Initial submission of first batch
                log.info(f"  Submitting initial batch (up to {max_pending_tasks} tasks)...")
                refill()

                log.info(f"  Initial batch submitted ({submitted_count} tasks). Processing with {n_jobs} workers...")
                log.info(f"  Note: Tasks have a {BASE_TASK_TIMEOUT // 60}-minute base timeout (doubled up to {MAX_RETRIES}x on retry). Per-future anti-deadlock active.")
                
                while pending_futures or retry_queue:
                    # Refill: after a batch of cancellations (e.g. all pending
                    # futures re-queued for retry) pending_futures may be empty
                    # while retry_queue still holds work to submit.
                    if not pending_futures and retry_queue:
                        refill()
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
                                
                                # Cancel all pending futures. No retry here: with
                                # the workers dead there is nothing to run them on,
                                # so every in-flight and not-yet-submitted task is
                                # recorded as a definitive failure (resumable via
                                # --resume failed_analyses.tsv).
                                for future in list(pending_futures):
                                    future.cancel()
                                    key, section_name, section_outdir = future_to_key[future]
                                    log.error(f"     ✗ CANCELLED (orphaned): {section_name}")
                                    record_failure(section_name, key, section_outdir, "orphaned", "worker death detected (likely OOM)")
                                    failed += 1
                                # Drain the retry queue and the main queue so the
                                # tasks that never ran are also logged for resume.
                                for task_tuple in list(retry_queue):
                                    _, _, _, _, name, key, _, _, _ = task_tuple
                                    record_failure(name, key, task_tuple[3], "orphaned", "worker death detected (likely OOM)")
                                    failed += 1
                                retry_queue.clear()
                                while True:
                                    try:
                                        task_tuple = next(task_iter)
                                    except StopIteration:
                                        break
                                    _, _, _, _, name, key, _, _, _ = task_tuple
                                    record_failure(name, key, task_tuple[3], "orphaned", "worker death detected (likely OOM)")
                                    failed += 1
                                submit_time.clear()
                                break
                            
                            last_worker_check = current_time
                        except Exception as e:
                            log.debug(f"Worker health check failed: {e}")
                    
                    # Per-future timeout: cancel ONLY the futures that exceeded
                    # their own budget (a single slow section must not cancel
                    # the whole batch). The SIGALRM in the worker should have
                    # killed the task at task_timeout; this is the backstop for
                    # futures that never report back (e.g. worker stuck in C
                    # code where signals don't fire).
                    now = time.time()
                    future_timeout = current_future_timeout()
                    for future in list(pending_futures):
                        t0 = submit_time.get(future)
                        if t0 is not None and (now - t0) > future_timeout:
                            future.cancel()
                            pending_futures.discard(future)
                            key, section_name, section_outdir = future_to_key[future]
                            submit_time.pop(future, None)
                            future_to_key.pop(future, None)
                            last_progress_time = time.time()
                            if handle_time_failure(key, section_name, section_outdir, "slow", f"exceeded {future_timeout}s"):
                                log.warning(f"  ↻ SLOW TASK ({completed+failed}/{len(tasks)}): {section_name} - exceeded {future_timeout}s, will retry")
                            else:
                                log.error(f"  ✗ SLOW TASK CANCELLED ({completed+failed}/{len(tasks)}): {section_name} - exceeded {future_timeout}s")

                    # Use short timeout on as_completed to check for stalls
                    if not pending_futures:
                        # All in-flight futures were just cancelled/re-queued;
                        # wait a beat before the next iteration resubmits them.
                        time.sleep(1)
                        continue
                    try:
                        done_iter = as_completed(pending_futures, timeout=5)
                        for future in done_iter:
                            pending_futures.discard(future)
                            submit_time.pop(future, None)
                            key, section_name, section_outdir = future_to_key[future]
                            future_to_key.pop(future, None)
                            
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
                                last_progress_time = time.time()
                                tto = task_timeout_by_key.get(key, BASE_TASK_TIMEOUT)
                                if handle_time_failure(key, section_name, section_outdir, "timeout", f"exceeded {tto}s"):
                                    log.warning(f"  ↻ TIMEOUT ({completed+failed}/{len(tasks)}): {section_name} - exceeded {tto}s, will retry")
                                else:
                                    log.error(f"  ✗ TIMEOUT ({completed+failed}/{len(tasks)}): {section_name} - exceeded {tto}s")
                            except Exception as e:
                                failed += 1
                                last_progress_time = time.time()
                                # Check if it's a worker crash (common patterns in error message)
                                if "process" in str(e).lower() and ("terminate" in str(e).lower() or "crash" in str(e).lower() or "abrupt" in str(e).lower()):
                                    log.error(f"  ✗ CRASH ({completed+failed}/{len(tasks)}): {section_name} - worker OOM or crash")
                                    record_failure(section_name, key, section_outdir, "crash", "worker OOM or crash")
                                else:
                                    log.error(f"  ✗ ERROR ({completed+failed}/{len(tasks)}): {section_name} - {type(e).__name__}")
                                    record_failure(section_name, key, section_outdir, "error", f"{type(e).__name__}: {e}")
                            
                            # CRITICAL: Submit next task to maintain max_pending_tasks in flight
                            # This lazy submission keeps memory usage constant regardless of total tasks.
                            # Retried tasks (retry_queue) are submitted before new ones.
                            if len(pending_futures) < max_pending_tasks:
                                try:
                                    submit_task(next_task())
                                    # Log every 50 submissions
                                    if submitted_count % 50 == 0:
                                        log.info(f"  → Submitted up to {submitted_count}/{len(tasks)} tasks (lazy mode, {len(retry_queue)} retries queued)")
                                except StopIteration:
                                    pass  # No more tasks to submit
                            
                            # Break inner loop to check stall timeout
                            break
                    except TimeoutError:
                        # No futures completed in 5 seconds, check for stall.
                        # This is the LAST-RESORT safety net (stall_timeout >
                        # future_timeout): it only fires when NO future has
                        # completed for longer than the per-future budget,
                        # i.e. a true deadlock (all workers stuck).
                        elapsed_since_progress = time.time() - last_progress_time
                        if elapsed_since_progress > current_stall_timeout():
                            log.error(f"  ✗ DEADLOCK DETECTED: No progress for {elapsed_since_progress:.0f}s")
                            log.error(f"     Cancelling {len(pending_futures)} remaining tasks to prevent infinite hang")
                            # Cancel all pending futures, retrying each with a
                            # doubled timeout (bounded by MAX_RETRIES).
                            any_retry = False
                            for future in list(pending_futures):
                                future.cancel()
                                pending_futures.discard(future)
                                key, section_name, section_outdir = future_to_key[future]
                                submit_time.pop(future, None)
                                future_to_key.pop(future, None)
                                if handle_time_failure(key, section_name, section_outdir, "deadlock", f"no progress for {elapsed_since_progress:.0f}s"):
                                    any_retry = True
                                else:
                                    log.error(f"     ✗ CANCELLED: {section_name}")
                            last_progress_time = time.time()
                            if not any_retry:
                                # Nothing left to retry: record the not-yet-run
                                # tasks (retry queue + main queue) and stop.
                                for task_tuple in list(retry_queue):
                                    _, _, _, _, name, key, _, _, _ = task_tuple
                                    record_failure(name, key, task_tuple[3], "deadlock", "no progress (stall)")
                                    failed += 1
                                retry_queue.clear()
                                while True:
                                    try:
                                        task_tuple = next(task_iter)
                                    except StopIteration:
                                        break
                                    _, _, _, _, name, key, _, _, _ = task_tuple
                                    record_failure(name, key, task_tuple[3], "deadlock", "no progress (stall)")
                                    failed += 1
                                submit_time.clear()
                                break
                
                if failed > 0:
                    log.warning(f"  {failed} tasks failed or cancelled out of {len(tasks)} total")
        else:
            # Sequential execution
            for df_source, cols, info, section_outdir, name, key, stat_test, bmk_filter, max_bmks in tasks:
                # Reconstruct DataFrame from pickled format
                if isinstance(df_source, tuple) and df_source[0] == 'pickled':
                    import pickle
                    df = pickle.loads(df_source[1])
                else:
                    df = pd.DataFrame(df_source)
                try:
                    results = analyze_section(df, cols, info, section_outdir, name, stat_test=stat_test, bmk_filter_cols=bmk_filter, max_bmks=max_bmks)
                except Exception as e:
                    log.error(f"  ✗ ERROR: {name} - {type(e).__name__}: {e}")
                    failed_analyses.append({
                        "vtype": vtype,
                        "mtype": key[0],
                        "section_key": key[1],
                        "section_name": name,
                        "section_outdir": section_outdir,
                        "reason": "error",
                        "detail": f"{type(e).__name__}: {e}",
                        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                    })
                    continue
                
                # Create slim results dict (same as parallel mode for consistency)
                # NOTE: must match the parallel-mode path (5_differential/), which is where
                # differential_analysis() actually writes differential_results.csv.
                slim_results = {
                    "differential_table": os.path.join(section_outdir, "5_differential", "differential_results.csv")
                }
                
                mtype, section_key = key
                all_results[vtype][mtype][section_key] = slim_results

        # ===============================================================
        # GLOBAL RANKING per mtype for this value_type
        # ===============================================================
        # Each mtype (aggregate / feature / sites) gets its own global_ranking dir:
        #   <vtype_dir>/<mtype>/global_ranking/
        # The second-pass analysis is run per mtype as well, so that the
        # significant BMKs are re-analysed within their own mtype context.
        for mtype in ("aggregate", "feature", "sites"):
            mdata = all_results[vtype].get(mtype, {})
            if not mdata:
                continue

            # The "feature-like" dataframe used to rebuild full rows in the
            # second pass: sites for the sites mtype, features otherwise.
            feat_like_df = sites_df if mtype == "sites" else feat_df

            gr_dir = os.path.join(vtype_dir, mtype, "global_ranking")
            log.info(f"\n--- GLOBAL RANKING for {vtype} / {mtype} ---")
            global_ranking({vtype: {mtype: mdata}}, gr_dir, stat_test=args.stat_test)

            # ===========================================================
            # 2nd Pass for only the top-ranked BMKs (per mtype)
            # ===========================================================
            # --- ANY comparison (union) ---
            df_sig_path = os.path.join(gr_dir, "significant_bmks_any_comparison", "significant_bmks.csv")
            if os.path.isfile(df_sig_path):
                df_sig = slurp_file(df_sig_path, separator=",")
                log.info(f"  Running second-pass analysis on any-comparison significant BMKs ({len(df_sig)} rows)")
                df_sig = filter_significant(df_sig, feat_like_df, agg_df, uid_col="uid", info=df_sig_path)
                output = os.path.join(gr_dir, "significant_bmks_any_comparison")
                results = analyze_section(df_sig, vcols, v_sample_info, output, "Section Significant", stat_test=args.stat_test, bmk_filter_cols=args.bmk_filter, max_bmks=args.max_bmks, enabled_tests=["descriptive", "multivariate", "correlation", "differential", "ranking", "classification", "stability", "heatmap"])

                # --- GLOBAL test only (omnibus) ---
                df_global_path = os.path.join(gr_dir, "significant_bmks_global_comparison", "significant_bmks.csv")
                if os.path.isfile(df_global_path):
                    df_global = slurp_file(df_global_path, separator=",")
                    log.info(f"  Running second-pass analysis on global-significant BMKs ({len(df_global)} rows)")
                    df_global = filter_significant(df_global, feat_like_df, agg_df, uid_col="uid", info=df_global_path)
                    output_global = os.path.join(gr_dir, "significant_bmks_global_comparison")
                    results = analyze_section(df_global, vcols, v_sample_info, output_global, "Section Significant", stat_test=args.stat_test, bmk_filter_cols=args.bmk_filter, max_bmks=args.max_bmks, enabled_tests=["descriptive", "multivariate", "correlation", "differential", "ranking", "classification", "stability", "heatmap"])

                # --- Per-pairwise comparison ---
                log.info(f"  Running second-pass analysis on per-condition significant BMKs from global ranking")
                base_dir = os.path.join(gr_dir, "significant_bmks_pairwise_comparison")
                all_sig_paths = glob.glob(os.path.join(base_dir, "**", "all_significant.csv"), recursive=True)
                log.info(f"Found {len(all_sig_paths)} files")
                for path in sorted(all_sig_paths):
                    folder = os.path.dirname(path)
                    folder_name = os.path.basename(folder)
                    df_sig = slurp_file(path, separator=",")
                    log.info(f"  Processing {folder_name} ({len(df_sig)} rows)")
                    df_sig = filter_significant(df_sig, feat_like_df, agg_df, uid_col="uid")
                    results = analyze_section(df_sig, vcols, v_sample_info, folder, "Section Significant", stat_test=args.stat_test, bmk_filter_cols=args.bmk_filter, max_bmks=args.max_bmks, enabled_tests=["descriptive", "multivariate", "correlation", "differential", "ranking", "classification", "stability", "heatmap"])
            else:
                log.warning(f"  No significant BMKs found in global ranking for {vtype} / {mtype}, skipping second-pass analysis")

            # ===========================================================
            # Save manifest per mtype (avoids overwriting across independent
            # runs; only the count relevant to this mtype is reported)
            # ===========================================================
            if mtype == "aggregate":
                counts = {"n_aggregates": len(agg_df) if agg_df is not None else 0}
            elif mtype == "feature":
                counts = {"n_features": len(feat_df) if feat_df is not None else 0}
            else:  # sites
                counts = {"n_sites": len(sites_df) if sites_df is not None else 0}

            manifest = {
                "value_types": value_types,
                "outdir": args.outdir,
                "mtype": mtype,
                **counts,
                "sample_info": sample_info,
            }
            with open(os.path.join(vtype_dir, mtype, "manifest.json"), "w") as f:
                json.dump(manifest, f, indent=2, default=str)

        # mtypes whose input file was provided but had no data rows: no
        # analysis ran, but the mtype dir + manifest are still written so
        # downstream consumers (Nextflow) find the expected output.
        for mtype in empty_mtypes:
            mtype_dir = os.path.join(vtype_dir, mtype)
            safe_mkdir(mtype_dir)
            if mtype == "aggregate":
                counts = {"n_aggregates": 0}
            elif mtype == "feature":
                counts = {"n_features": 0}
            else:  # sites
                counts = {"n_sites": 0}
            manifest = {
                "value_types": value_types,
                "outdir": args.outdir,
                "mtype": mtype,
                **counts,
                "sample_info": sample_info,
            }
            with open(os.path.join(mtype_dir, "manifest.json"), "w") as f:
                json.dump(manifest, f, indent=2, default=str)
            log.info(f"  Wrote empty mtype dir + manifest: {mtype_dir}")

    # Write the failed/cancelled analyses log so the user can re-run them later.
    if failed_analyses:
        failed_log_path = os.path.join(args.outdir, "failed_analyses.tsv")
        pd.DataFrame(failed_analyses).to_csv(failed_log_path, sep="\t", index=False)
        log.warning(f"\n{len(failed_analyses)} analyses failed or were cancelled. See {failed_log_path} to re-run them.")
    else:
        # Remove a stale log from a previous run so it doesn't mislead.
        stale = os.path.join(args.outdir, "failed_analyses.tsv")
        if os.path.isfile(stale):
            os.remove(stale)

    log.info(f"\nAnalysis complete. Results saved to {args.outdir}/")


if __name__ == "__main__":
    main()
