"""barometer.utils – Data loading, UID creation, sample-column parsing, helpers.

Shared low-level utilities used by analysis.py, ranking.py and __main__.py.
No dependency on the other barometer modules (leaf of the import graph).
"""

import logging
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# DataFrame helpers
# ---------------------------------------------------------------------------

# to avoid DtypeWarning
def slurp_file(path, sample_cols=None, separator=None):
    if sample_cols:
        feat_df = pd.read_csv(
            path,
            sep=separator,
            dtype=str,
            low_memory=False
        )
        # convertir uniquement les colonnes numériques
        feat_df[sample_cols] = feat_df[sample_cols].apply(
            pd.to_numeric, errors="coerce"
        )
    else:
        feat_df = pd.read_csv(
            path,
            sep=separator,
            dtype={"SeqID": str, "Start": str, "End": str, "Strand": str},
            low_memory=False
        )

    return feat_df


def filter_significant(df_sig, feat_df, agg_df, uid_col="uid", info=None):
    """Filter df to only include rows where uid is in sig_df[uid_col]."""
    # ── get significant ─────────────────────────────────────────────
    if uid_col not in df_sig.columns:
        log.error(f"Column '{uid_col}' not found in DataFrame: {info}")
        return None
    else:
        sig_uids = df_sig[uid_col].tolist()

    # ── Filtrer les significatifs ─────────────────────────────────────────────
    feat_df_sig = feat_df[feat_df[uid_col].isin(sig_uids)] if feat_df is not None and len(feat_df) > 0 else None
    agg_df_sig  = agg_df[agg_df[uid_col].isin(sig_uids)]   if agg_df  is not None and len(agg_df)  > 0 else None

    # ── Jointure  ─────────────────────────────────────────────────
    if feat_df_sig is not None and agg_df_sig is not None:
        # Les deux existent → on fusionne
        common_cols = set(agg_df_sig.columns) & set(feat_df_sig.columns)
        # garder seulement 'uid' comme clé commune
        cols_to_drop = [c for c in common_cols if c != uid_col]
        feat_df_sig_clean = feat_df_sig.drop(columns=cols_to_drop)
        df_sig = pd.merge(agg_df_sig, feat_df_sig_clean, on=uid_col, how="outer")
        log.info(f"  Merged feat + agg → {len(df_sig)} rows")

    elif feat_df_sig is not None:
        # Seulement feat
        df_sig = feat_df_sig
        log.info(f"  Only feat_df available → {len(df_sig)} rows")

    elif agg_df_sig is not None:
        # Seulement agg
        df_sig = agg_df_sig
        log.info(f"  Only agg_df available → {len(df_sig)} rows")

    else:
        # Aucun des deux
        log.warning("  No data available (feat_df and agg_df are both empty/None)")
        df_sig = pd.DataFrame()

    #print all column
    log.debug(f"  Final df_sig columns: {df_sig.columns.tolist()}")

    return df_sig


def safe_mkdir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def harmonize_columns(df):
    """Replace '..' par '::' dfor compatibility in case it has been proceceed by R."""
    df.columns = [c.replace('..', '::') for c in df.columns]
    return df


def prepare_df_for_task(df):
    """Convert DataFrame to a memory-efficient format for task submission.

    Uses pickle which is much faster than .to_dict('list') for large DataFrames.
    Returns a tuple ('pickled', bytes_data) that can be unpickled in the worker.
    """
    import pickle
    return ('pickled', pickle.dumps(df, protocol=pickle.HIGHEST_PROTOCOL))


# work line by line (Series)
def create_uid_column(df=None, meta_cols=["SeqID", "ParentIDs", "ID", "Mtype", "Ptype", "Type", "Ctype", "Mode"]):
    """Create a unique identifier (vectorized, robust)"""

    def clean_part(series):
        cleaned = series.astype(str).str.strip()
        # remove leading ., or ., e.g. for ParentsID
        cleaned = cleaned.str.replace(r"^\.?,?", "", regex=True)
        cleaned = cleaned.replace(r"^\.*$", pd.NA, regex=True)
        cleaned = cleaned.replace(r"^\s*$", pd.NA, regex=True)
        return cleaned.fillna("").astype(str)

    # ── Mtype (prefix) ───────────────────────────────────
    uid = pd.Series("", index=df.index)

    if "Mtype" in df.columns:
        mtype_clean = df["Mtype"].astype(str).str.lower().str.strip()

        # Utiliser loc avec le masque booléen
        uid.loc[mtype_clean == "aggregate"] = "agg"
        uid.loc[mtype_clean == "feature"] = "feat"

    composite_parts = []

    # ── Type ─────────────────────────────────────────────
    if "Type" in df.columns:
        composite_parts.append(clean_part(df["Type"]))

    # ── SeqID ────────────────────────────────────────────
    if "SeqID" in df.columns:
        composite_parts.append(clean_part(df["SeqID"]))

    # ── Ptype ────────────────────────────────────────────
    if "Ptype" in df.columns:
        composite_parts.append(clean_part(df["Ptype"]))

    # ── ParentIDs ────────────────────────────────────────
    if "ParentIDs" in df.columns:
        composite_parts.append(clean_part(df["ParentIDs"]))

    # ── ID → seulement si Type n’est pas sequence/global ──
    if "ID" in df.columns and "Type" in df.columns:
        type_clean = df["Type"].astype(str).str.lower().str.strip()
         # garder ID si Type NOT IN ["global", "sequence"]
        id_allowed = df["ID"].where(~type_clean.isin(["global", "sequence"]), pd.NA)
        composite_parts.append(clean_part(id_allowed))

    # ── Ctype ───────────────────────────────────────────
    if "Ctype" in df.columns:
        composite_parts.append(clean_part(df["Ctype"]))

    # ── Mode ────────────────────────────────────────────
    if "Mode" in df.columns:
        composite_parts.append(clean_part(df["Mode"]))

    # ── Concaténation des parties
    for part in composite_parts:
        uid = uid + "-" + part

    # Supprimer doubles tirets et tirets devant/derrière
    uid = uid.str.replace(r"-+", "-", regex=True).str.strip("-")

    # Mettre uid en première colonne
    df.insert(0, "uid", uid)  # in-place, no copy needed (df is a fresh DataFrame from slurp_file)

    return df


# ---------------------------------------------------------------------------
# Sample-column parsing
# ---------------------------------------------------------------------------

def parse_sample_columns(columns):
    """Parse sample column names like 'test1::rain_chr21_small::rep1::espf'.

    Returns a list of dicts with keys: col, group, sample, rep, value_type.
    """
    parsed = []
    meta_cols = {"uid", "SeqID", "ParentIDs", "ID", "Mtype", "Ptype", "Type", "Ctype", "Mode", "Start", "End", "Strand"}
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


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------

def numeric_df(df, cols):
    """Return a copy with the given columns cast to numeric (coerce errors)."""
    out = df.copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def save_fig(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
