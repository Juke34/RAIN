#!/usr/bin/env python3
"""
barometer_report.py – Generate an interactive HTML report from barometer_analyze.py results.

Assembles all tables, figures, and analyses into a multi-tab HTML report using
Jinja2 templating with interactive DataTables, Plotly for image zoom, and
collapsible sections for navigation.

Usage:
    python barometer_report.py [--results DIR] [--output FILE]
"""

import argparse
import base64
import csv
import glob
import json
import logging
import os
import re
import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd
from jinja2 import Template

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EMBED_IMAGES = False
RESULTS_DIR = "llmeter_results"


def img_src(path):
    """Return an image src: base64 data URI if embedding, else relative path."""
    if not os.path.exists(path):
        return ""
    if EMBED_IMAGES:
        with open(path, "rb") as f:
            data = f.read()
        ext = path.rsplit(".", 1)[-1].lower()
        mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
                "svg": "image/svg+xml"}.get(ext, "image/png")
        return f"data:{mime};base64,{base64.b64encode(data).decode()}"
    else:
        # Return relative path from report file location
        return path


def csv_to_html_table(path, max_rows=500, table_id=None):
    """Read a CSV and return an HTML table string."""
    if not os.path.exists(path):
        return ""
    try:
        df = pd.read_csv(path)
    except Exception:
        return ""
    if len(df) > max_rows:
        df = df.head(max_rows)
    tid = f' id="{table_id}"' if table_id else ""
    classes = "display compact stripe hover"
    html = df.to_html(index=False, classes=classes, border=0, na_rep="NA",
                      float_format=lambda x: f"{x:.4g}" if abs(x) > 1e-4 else f"{x:.2e}")
    html = html.replace("<table ", f"<table{tid} ", 1) if table_id else html
    return html


def find_images(directory):
    """Find all PNG images in a directory."""
    imgs = []
    for ext in ["*.png", "*.jpg", "*.jpeg"]:
        imgs.extend(sorted(glob.glob(os.path.join(directory, ext))))
    return imgs


def find_csvs(directory):
    """Find all CSV files in a directory."""
    return sorted(glob.glob(os.path.join(directory, "*.csv")))


def sanitize_id(text):
    """Create a safe HTML id from text."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", str(text))


def prettify_section_name(rel_path):
    """Convert a path like '1_global/1_2_by_ctype_all_isoforms' to 'Global / All Isoforms'."""
    parts = rel_path.split(os.sep)
    clean_parts = []
    for part in parts:
        # Remove leading numeric prefix like "1_", "1_1_", "2_3_", etc.
        part = re.sub(r'^\d+(_\d+)*_', '', part)
        # Remove structural grouping prefixes
        part = re.sub(r'^by_[a-z]+_', '', part)
        # Remove leading underscore (e.g., _all → all)
        part = re.sub(r'^_', '', part)
        # Replace underscores with spaces
        part = part.replace("_", " ").strip()
        # Remove trailing " Agg" or " agg" suffixes (e.g. "feature_agg" → "Feature")
        part = re.sub(r'\s+[Aa]gg$', '', part)
        # Capitalise the first letter of each word without lowercasing the rest
        # (preserves biological abbreviations like ncRNA, snRNA, etc.)
        part = re.sub(r'\b\w', lambda m: m.group().upper(), part)
        if part:
            clean_parts.append(part)
    # Merge pure-numeric parts into the previous component ("Chr" + "21" → "Chr 21")
    merged = []
    for p in clean_parts:
        if p.isdigit() and merged:
            merged[-1] = merged[-1] + " " + p
        else:
            merged.append(p)
    result = " / ".join(merged)
    # Strip leading "All " — these sections live inside an "All" subgroup already
    if result.startswith("All "):
        result = result[4:]
    return result


# Group-name mapping for first-level directory components
_GROUP_MAP = {
    "global":  "Global",
    "sequence": "Sequence",
    "feature": "Feature",
}


def prettify_group_name(dirname):
    """Return a human-readable group name for the first directory component."""
    if dirname in _GROUP_MAP:
        return _GROUP_MAP[dirname]
    # For feature-mtype top-level types (gene, ncRNA_gene, …)
    name = dirname.replace("_", " ").strip()
    name = re.sub(r'\s+[Aa]gg$', '', name)
    name = re.sub(r'\b\w', lambda m: m.group().upper(), name)
    return name


_CTYPE_RE = re.compile(
    r'(?:^|_)by_ptype_(.+?)_(all_isoforms|chimaera|longest_isoform)$'
)

_MODE_DISPLAY = {
    "all_isoforms":      "Isoforms",
    "chimaera":          "Chimaera",
    "longest_isoform":   "Longest Isoform",
    "feature":           "Feature",
}


def _read_section_meta(section_dir):
    """
    Read Ptype, Ctype, Mode from the first data row of data.csv in section_dir.
    Returns a dict with keys 'Ptype', 'Ctype', 'Mode', or {} if not found.
    """
    csv_path = os.path.join(section_dir, "data.csv")
    if not os.path.exists(csv_path):
        return {}
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            row = next(reader, None)
            if row and "Ptype" in row and "Ctype" in row and "Mode" in row:
                return {"Ptype": row["Ptype"], "Ctype": row["Ctype"], "Mode": row["Mode"]}
    except Exception:
        pass
    return {}


def _prettify_mode(mode_raw):
    """Human-readable mode label with fallback."""
    if mode_raw in _MODE_DISPLAY:
        return _MODE_DISPLAY[mode_raw]
    name = mode_raw.replace("_", " ")
    return re.sub(r'\b\w', lambda m: m.group().upper(), name)


def _extract_ptype(dirname):
    """Return (ptype_raw, ctype_raw) from a by_ptype dirname, or (None, None)."""
    m = _CTYPE_RE.search(dirname)
    if m:
        return m.group(1), m.group(2)
    return None, None


def _prettify_ptype(ptype_raw):
    """Human-readable ptype label: 'ncRNA_gene' → 'NcRNA Gene'."""
    name = ptype_raw.replace("_", " ")
    name = re.sub(r'\b\w', lambda m: m.group().upper(), name)
    return name


# ---------------------------------------------------------------------------
# Build report structure from results directory
# ---------------------------------------------------------------------------

def build_section_content(section_dir, section_name, table_counter):
    """Build HTML content for a single analysis section."""
    content = ""
    subdirs = ["qc", "descriptive", "differential", "multivariate",
               "correlation", "ranking", "classification", "stability", "batch"]

    analysis_labels = {
        "qc": "1. Quality Control",
        "descriptive": "3. Descriptive Statistics",
        "differential": "4. Differential Analysis & 5. Multiple Testing Correction",
        "multivariate": "6. Multivariate Analysis (PCA)",
        "correlation": "7. Correlation / Network Analysis",
        "ranking": "8. Feature Selection / Biomarker Ranking",
        "classification": "9. Classification / Predictive Modeling",
        "stability": "12. Stability / Robustness Analysis",
        "batch": "11. Batch Effect Detection",
    }

    # Section-level heatmap
    heatmap_path = os.path.join(section_dir, "heatmap.png")
    if os.path.exists(heatmap_path):
        b64 = img_src(heatmap_path)
        content += f'<div class="figure-container"><h5>13. Visualization – Overview Heatmap</h5>'
        content += f'<img src="{b64}" class="report-img zoomable" alt="Heatmap"></div>\n'

    for subdir in subdirs:
        sub_path = os.path.join(section_dir, subdir)
        if not os.path.isdir(sub_path):
            continue

        label = analysis_labels.get(subdir, subdir.title())
        content += f'<div class="analysis-block"><h5>{label}</h5>\n'

        # Tables
        for csv_path in find_csvs(sub_path):
            fname = os.path.basename(csv_path).replace(".csv", "")
            table_counter[0] += 1
            tid = f"dt_{table_counter[0]}"
            content += f'<div class="table-section"><h6>{fname.replace("_", " ").title()}</h6>\n'
            content += f'<div class="table-responsive">{csv_to_html_table(csv_path, table_id=tid)}</div>\n'
            content += f'<script>$(document).ready(function(){{ if($.fn.DataTable){{ try{{ $("#{tid}").DataTable({{paging:true,pageLength:20,searching:true,ordering:true,scrollX:true}});}}catch(e){{}}}} }});</script>\n'
            content += '</div>\n'

        # Images
        for img_path in find_images(sub_path):
            fname = os.path.basename(img_path).replace(".png", "").replace("_", " ").title()
            b64 = img_src(img_path)
            if b64:
                content += f'<div class="figure-container"><h6>{fname}</h6>'
                content += f'<img src="{b64}" class="report-img zoomable" alt="{fname}"></div>\n'

        content += '</div>\n'

    return content


def build_report_data(results_dir):
    """Walk the results directory and build the report structure."""
    manifest_path = os.path.join(results_dir, "manifest.json")
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)

    value_types = manifest.get("value_types", [])
    if not value_types:
        # Discover from directory
        for d in sorted(os.listdir(results_dir)):
            dp = os.path.join(results_dir, d)
            if os.path.isdir(dp) and d not in ("global_ranking",):
                value_types.append(d)

    report = OrderedDict()
    table_counter = [0]

    for vtype in value_types:
        vtype_dir = os.path.join(results_dir, vtype)
        if not os.path.isdir(vtype_dir):
            continue

        report[vtype] = OrderedDict()

        # Aggregate and Feature tabs
        for mtype in ["aggregate", "feature"]:
            mtype_dir = os.path.join(vtype_dir, mtype)
            if not os.path.isdir(mtype_dir):
                continue

            report[vtype][mtype] = OrderedDict()
            # Walk all section directories
            for root, dirs, files in os.walk(mtype_dir):
                # Only process leaf directories that contain actual results
                if any(f.endswith(".csv") or f.endswith(".png") for f in files):
                    rel = os.path.relpath(root, mtype_dir)
                    parts = rel.split(os.sep)
                    first = parts[0]
                    group_name = prettify_group_name(first)

                    # Detect subgroup / sub-subgroup
                    ssg_key  = None
                    ssg_name = None
                    if first == "global" and len(parts) >= 2:
                        sec_dir   = parts[1]
                        remainder = os.sep.join(parts[1:])
                        if "by_ptype_" not in sec_dir:
                            # all_sites / by_ctype_* → "All" subgroup
                            sg_key  = "all"
                            sg_name = "All"
                            section_name = prettify_section_name(remainder)
                        else:
                            # by_ptype_<ptype>_<mode> → subgroup per ptype
                            ptype_raw, ctype_raw = _extract_ptype(sec_dir)
                            if ptype_raw:
                                sg_key       = f"ptype_{ptype_raw}"
                                sg_name      = _prettify_ptype(ptype_raw)
                                section_name = _prettify_mode(ctype_raw) if ctype_raw else prettify_section_name(remainder)
                            else:
                                sg_key       = ""
                                sg_name      = ""
                                section_name = prettify_section_name(remainder)
                    elif first == "sequence" and len(parts) >= 3:
                        sg_raw  = parts[1]          # e.g. "21" or "all_sequence_bmks"
                        sg_name = "All Sequences" if sg_raw == "all_sequence_bmks" else f"Chr {sg_raw}"
                        sg_key  = sg_raw
                        sec_dir   = parts[2]
                        remainder = os.sep.join(parts[2:])
                        if "by_ptype_" not in sec_dir:
                            # all_sites / by_ctype_* → "All" sub-subgroup
                            ssg_key  = "all"
                            ssg_name = "All"
                            section_name = prettify_section_name(remainder)
                        else:
                            # by_ptype_<ptype>_<mode> → sub-subgroup per ptype
                            ptype_raw, ctype_raw = _extract_ptype(sec_dir)
                            if ptype_raw:
                                ssg_key      = f"ptype_{ptype_raw}"
                                ssg_name     = _prettify_ptype(ptype_raw)
                                section_name = _prettify_mode(ctype_raw) if ctype_raw else prettify_section_name(remainder)
                            else:
                                ssg_key  = None
                                ssg_name = None
                                section_name = prettify_section_name(remainder)
                    elif first == "feature" and len(parts) >= 3 and parts[1] == "all_feature_together":
                        # all_feature_together/{safe_name} — all chromosomes pooled
                        sec_dir     = parts[2]
                        section_dir = os.path.join(mtype_dir, first, parts[1], sec_dir)
                        meta        = _read_section_meta(section_dir)
                        if meta:
                            ptype_raw    = meta["Ptype"]
                            ctype_raw    = meta["Ctype"]
                            mode         = meta["Mode"]
                            sg_key       = ptype_raw
                            sg_name      = _prettify_ptype(ptype_raw)
                            ssg_key      = ctype_raw
                            ssg_name     = ctype_raw.replace("_", " ").replace("-", " ").title()
                            section_name = _prettify_mode(mode)
                        else:
                            sg_key       = "all_feature_together"
                            sg_name      = "All Features Together"
                            ssg_key      = None
                            ssg_name     = None
                            section_name = prettify_section_name(parts[2])
                    elif first == "feature" and len(parts) >= 4 and parts[1] == "by_sequence":
                        # by_sequence/{chrom}/{safe_name} — features grouped by chromosome
                        chrom       = parts[2]
                        sec_dir     = parts[3]
                        section_dir = os.path.join(mtype_dir, first, parts[1], chrom, sec_dir)
                        meta        = _read_section_meta(section_dir)
                        sg_key      = f"seq_{chrom}"
                        sg_name     = f"Chr {chrom}"
                        if meta:
                            ptype_raw    = meta["Ptype"]
                            ctype_raw    = meta["Ctype"]
                            mode         = meta["Mode"]
                            ssg_key      = ctype_raw
                            ssg_name     = ctype_raw.replace("_", " ").replace("-", " ").title()
                            section_name = f"{_prettify_ptype(ptype_raw)} – {_prettify_mode(mode)}"
                        else:
                            ssg_key      = None
                            ssg_name     = None
                            section_name = prettify_section_name(parts[3])
                    elif first == "feature" and len(parts) >= 2:
                        # Flat fallback (old structure)
                        sec_dir      = parts[1]
                        section_dir  = os.path.join(mtype_dir, first, sec_dir)
                        meta         = _read_section_meta(section_dir)
                        if meta:
                            ptype_raw    = meta["Ptype"]
                            ctype_raw    = meta["Ctype"]
                            mode         = meta["Mode"]
                            sg_key       = ptype_raw
                            sg_name      = _prettify_ptype(ptype_raw)
                            ssg_key      = ctype_raw
                            ssg_name     = ctype_raw.replace("_", " ").replace("-", " ").title()
                            section_name = _prettify_mode(mode)
                        else:
                            sg_key       = ""
                            sg_name      = ""
                            section_name = prettify_section_name(os.sep.join(parts[1:]))
                    else:
                        sg_key   = ""
                        sg_name  = ""
                        remainder = os.sep.join(parts[1:]) if len(parts) > 1 else ""
                        section_name = prettify_section_name(remainder) if remainder else prettify_section_name(first)

                    # Ensure group exists
                    if group_name not in report[vtype][mtype]:
                        report[vtype][mtype][group_name] = OrderedDict()

                    # Ensure subgroup bucket exists
                    # Structure: report[vtype][mtype][group_name][sg_key] = {
                    #   "name": sg_name, "sections": OrderedDict(), "subgroups": OrderedDict()}
                    grp = report[vtype][mtype][group_name]
                    if sg_key not in grp:
                        grp[sg_key] = {"name": sg_name, "sections": OrderedDict(), "subgroups": OrderedDict()}

                    content = build_section_content(root, section_name, table_counter)
                    if content.strip():
                        if (first in ("sequence", "feature")) and ssg_key:
                            # Store in a sub-subgroup (All / Gene / …) inside the chr subgroup
                            ssg = grp[sg_key]["subgroups"]
                            if ssg_key not in ssg:
                                ssg[ssg_key] = {"name": ssg_name, "sections": OrderedDict()}
                            ssg[ssg_key]["sections"][rel] = {"name": section_name, "content": content}
                        else:
                            grp[sg_key]["sections"][rel] = {"name": section_name, "content": content}

    # Global ranking
    gr_dir = os.path.join(results_dir, "global_ranking")
    global_content = ""
    if os.path.isdir(gr_dir):
        for csv_path in find_csvs(gr_dir):
            fname = os.path.basename(csv_path).replace(".csv", "")
            table_counter[0] += 1
            tid = f"dt_{table_counter[0]}"
            global_content += f'<div class="table-section"><h5>{fname.replace("_", " ").title()}</h5>\n'
            global_content += f'<div class="table-responsive">{csv_to_html_table(csv_path, table_id=tid)}</div>\n'
            global_content += f'<script>$(document).ready(function(){{ if($.fn.DataTable){{ try{{ $("#{tid}").DataTable({{paging:true,pageLength:30,searching:true,ordering:true,scrollX:true}});}}catch(e){{}}}} }});</script>\n'
            global_content += '</div>\n'

        for img_path in find_images(gr_dir):
            fname = os.path.basename(img_path).replace(".png", "").replace("_", " ").title()
            b64 = img_src(img_path)
            if b64:
                global_content += f'<div class="figure-container"><h6>{fname}</h6>'
                global_content += f'<img src="{b64}" class="report-img zoomable" alt="{fname}"></div>\n'

    return report, manifest, global_content


# ---------------------------------------------------------------------------
# HTML Template
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>barometer Biomarker Analysis Report</title>
<!-- jQuery -->
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<!-- DataTables -->
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<!-- Bootstrap 5 -->
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js"></script>
<style>
:root {
    --sidebar-width: 280px;
}
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    margin: 0;
    background: #f8f9fa;
}
/* Sidebar */
#sidebar {
    position: fixed;
    top: 0;
    left: 0;
    width: var(--sidebar-width);
    height: 100vh;
    background: #2c3e50;
    color: white;
    overflow-y: auto;
    padding-top: 60px;
    z-index: 1000;
    transition: transform 0.3s;
}
#sidebar .nav-header {
    padding: 15px 20px;
    font-size: 1.1em;
    font-weight: bold;
    color: #ecf0f1;
    border-bottom: 1px solid #34495e;
}
#sidebar .nav-link {
    color: #bdc3c7;
    padding: 6px 20px;
    font-size: 0.85em;
    cursor: pointer;
    display: block;
    text-decoration: none;
}
#sidebar .nav-link:hover, #sidebar .nav-link.active {
    color: #fff;
    background: #34495e;
}
#sidebar .nav-link.level-1 { padding-left: 20px; font-weight: bold; color: #3498db; }
#sidebar .nav-link.level-2 { padding-left: 35px; font-size: 0.82em; }
#sidebar .nav-link.level-2g { padding-left: 44px; font-size: 0.82em; color: #abebc6; font-style: italic; }
#sidebar .nav-link.level-3 { padding-left: 56px; font-size: 0.76em; color: #85c1e9; font-style: italic; }
#sidebar .nav-link.level-4 { padding-left: 70px; font-size: 0.72em; }
#sidebar .nav-link.level-5 { padding-left: 84px; font-size: 0.68em; color: #a9dfbf; }
/* Hide deep levels in sidebar by default - aggregate: hide >=4, feature: hide >=3 */
#sidebar .nav-link.sidebar-collapsed {
    display: none;
}
/* Collapsible sidebar items with hidden children */
#sidebar .nav-link.has-hidden-children::before {
    content: "▸ ";
    font-family: monospace;
    margin-right: 3px;
}
#sidebar .nav-link.has-hidden-children.expanded::before {
    content: "▾ ";
}
/* Main content */
#main-content {
    margin-left: var(--sidebar-width);
    padding: 20px 30px;
    padding-top: 70px;
}
/* Top bar */
#topbar {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    height: 55px;
    background: #1a252f;
    color: white;
    display: flex;
    align-items: center;
    padding: 0 20px 0 calc(var(--sidebar-width) + 20px);
    z-index: 1001;
    font-size: 1.2em;
    font-weight: bold;
    box-shadow: 0 2px 5px rgba(0,0,0,0.2);
}
/* Tabs */
.nav-tabs .nav-link { font-size: 0.9em; }
.tab-content { padding-top: 15px; }
/* Tables */
.table-responsive { margin-bottom: 15px; overflow-x: auto; }
table.dataTable { font-size: 0.8em; width: 100% !important; }
table.dataTable th, table.dataTable td { white-space: nowrap; padding: 4px 8px !important; }
.table-section { margin-bottom: 20px; }
/* Figures */
.figure-container { margin: 10px 0 20px 0; text-align: center; }
.report-img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }
.report-img.zoomable { cursor: zoom-in; transition: transform 0.3s; }
/* Analysis blocks */
.analysis-block {
    border-left: 3px solid #3498db;
    padding-left: 15px;
    margin-bottom: 20px;
}
.analysis-block h5 { color: #2c3e50; margin-top: 10px; }
/* Section cards */
.section-card {
    background: white;
    border-radius: 8px;
    padding: 20px;
    margin-bottom: 20px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}
.section-card h3 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 8px; }
.section-card h4 { color: #34495e; margin-top: 15px; }
/* Section collapse */
.section-toggle {
    cursor: pointer;
    user-select: none;
}
.section-toggle::before {
    content: "▸ ";
    font-family: monospace;
}
.section-toggle.open::before {
    content: "▾ ";
}
/* Image zoom modal */
#img-modal {
    display: none;
    position: fixed;
    top: 0; left: 0; right: 0; bottom: 0;
    background: rgba(0,0,0,0.85);
    z-index: 2000;
    cursor: zoom-out;
    padding: 40px;
    text-align: center;
}
#img-modal img {
    max-width: 95%;
    max-height: 90vh;
    border-radius: 4px;
}
/* Mode filter */
.mode-filter { margin-bottom: 10px; }
.mode-filter label { margin-right: 15px; font-size: 0.85em; }
/* Group blocks */
.group-block {
    background: #eaf2f8;
    border-radius: 8px;
    padding: 10px 20px 5px;
    margin-bottom: 15px;
    border-left: 4px solid #2980b9;
}
.group-toggle {
    cursor: pointer;
    user-select: none;
    color: #1a5276;
    margin: 0;
    padding: 6px 0;
    font-size: 1.15em;
    font-weight: 600;
}
.group-toggle::before { content: "▸ "; font-family: monospace; }
.group-toggle.open::before { content: "▾ "; }
.group-body { margin-top: 10px; }
/* Subgroup blocks (e.g. Chr 21 inside Sequence) */
.subgroup-block {
    background: #d6eaf8;
    border-radius: 6px;
    padding: 8px 15px 4px;
    margin-bottom: 10px;
    border-left: 3px solid #5dade2;
}
.subgroup-toggle {
    cursor: pointer;
    user-select: none;
    color: #1f618d;
    margin: 0;
    padding: 5px 0;
    font-size: 1em;
    font-weight: 600;
}
.subgroup-toggle::before { content: "▸ "; font-family: monospace; }
.subgroup-toggle.open::before { content: "▾ "; }
.subgroup-body { margin-top: 8px; }
/* Sub-subgroup blocks (e.g. "All" inside Chr 21) */
.subsubgroup-block {
    background: #d4efdf;
    border-radius: 5px;
    padding: 6px 12px 3px;
    margin-bottom: 8px;
    border-left: 3px solid #52be80;
}
.subsubgroup-toggle {
    cursor: pointer;
    user-select: none;
    color: #1e8449;
    margin: 0;
    padding: 4px 0;
    font-size: 0.9em;
    font-weight: 600;
}
.subsubgroup-toggle::before { content: "▸ "; font-family: monospace; }
.subsubgroup-toggle.open::before { content: "▾ "; }
.subsubgroup-body { margin-top: 6px; }
/* Print */
@media print {
    #sidebar, #topbar { display: none !important; }
    #main-content { margin-left: 0; padding-top: 0; }
}
</style>
</head>
<body>

<div id="topbar">barometer Biomarker Analysis Report</div>

<nav id="sidebar">
    <div class="nav-header">Navigation</div>
    <a class="nav-link level-1" href="#global-ranking">🏆 Global Ranking</a>
    {% for vtype, mtype_data in report.items() %}
    <a class="nav-link level-1" href="#vtype-{{ sanitize(vtype) }}">{{ vtype }}</a>
    {% for mtype, groups in mtype_data.items() %}
    <a class="nav-link level-2" href="#mtype-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}" data-mtype="{{ mtype }}">{{ mtype|title }}</a>
    {% for group_name, subgroups in groups.items() %}
    <a class="nav-link level-2g" href="#group-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(group_name) }}" data-mtype="{{ mtype }}">{{ group_name }}</a>
    {% for sg_key, sg_data in subgroups.items() %}
      {% if sg_key != "" %}
    <a class="nav-link level-3" href="#sg-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(group_name) }}-{{ sanitize(sg_key) }}" data-mtype="{{ mtype }}">{{ sg_data.name }}</a>
      {% endif %}
      {% for ssg_key, ssg_data in sg_data.subgroups.items() %}
    <a class="nav-link {% if sg_key != '' %}level-4{% else %}level-3{% endif %}{% if mtype == 'aggregate' and sg_key != '' %} sidebar-collapsed{% elif mtype == 'feature' %} sidebar-collapsed{% endif %}" href="#ssg-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sg_key) }}-{{ sanitize(ssg_key) }}" data-mtype="{{ mtype }}">{{ ssg_data.name }}</a>
        {% for sec_key, sec_data in ssg_data.sections.items() %}
    <a class="nav-link {% if sg_key != '' %}level-5{% else %}level-4{% endif %}{% if mtype == 'aggregate' %} sidebar-collapsed{% elif mtype == 'feature' %} sidebar-collapsed{% endif %}" href="#sec-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sec_key) }}" data-mtype="{{ mtype }}">{{ sec_data.name }}</a>
        {% endfor %}
      {% endfor %}
      {% for sec_key, sec_data in sg_data.sections.items() %}
    <a class="nav-link {% if sg_key != '' %}level-4{% else %}level-3{% endif %}{% if mtype == 'aggregate' and sg_key != '' %} sidebar-collapsed{% elif mtype == 'feature' %} sidebar-collapsed{% endif %}" href="#sec-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sec_key) }}" data-mtype="{{ mtype }}">{{ sec_data.name }}</a>
      {% endfor %}
    {% endfor %}
    {% endfor %}
    {% endfor %}
    {% endfor %}
</nav>

<div id="main-content">

    <!-- Summary -->
    <div class="section-card">
        <h3>Report Summary</h3>
        <p><strong>Aggregates:</strong> {{ manifest.get('n_aggregates', 'N/A') }} rows &nbsp;|&nbsp;
           <strong>Features:</strong> {{ manifest.get('n_features', 'N/A') }} rows &nbsp;|&nbsp;
           <strong>Value types:</strong> {{ manifest.get('value_types', [])|join(', ') }}</p>
        <p><strong>Samples:</strong></p>
        <ul>
        {% for s in manifest.get('sample_info', [])[:20] %}
            <li>{{ s.col }} (group: {{ s.group }}, sample: {{ s.sample }}, rep: {{ s.rep }})</li>
        {% endfor %}
        </ul>
        <p><em>Analyses performed: QC, Descriptive Statistics, Differential Editing (Mann-Whitney U, Kruskal-Wallis, ANOVA),
        Multiple Testing Correction (FDR-BH), PCA, Hierarchical Clustering, Correlation Analysis, Random Forest Feature
        Importance, Classification (RF, GBT, LDA), Stability/Robustness (CV, replicate concordance), Batch Effect Detection.</em></p>
    </div>

    <!-- Global Ranking -->
    <div class="section-card" id="global-ranking">
        <h3>🏆 Global Biomarker Ranking</h3>
        {{ global_content|safe }}
    </div>

    <!-- Value type tabs -->
    <ul class="nav nav-tabs" id="vtypeTabs" role="tablist">
    {% for vtype in report.keys() %}
        <li class="nav-item">
            <button class="nav-link {% if loop.first %}active{% endif %}"
                    id="tab-{{ sanitize(vtype) }}" data-bs-toggle="tab"
                    data-bs-target="#pane-{{ sanitize(vtype) }}" type="button">{{ vtype }}</button>
        </li>
    {% endfor %}
    </ul>

    <div class="tab-content" id="vtypeContent">
    {% for vtype, mtype_data in report.items() %}
        <div class="tab-pane fade {% if loop.first %}show active{% endif %}"
             id="pane-{{ sanitize(vtype) }}" role="tabpanel">

            <div id="vtype-{{ sanitize(vtype) }}"></div>

            <!-- MType tabs -->
            <ul class="nav nav-tabs" id="mtypeTabs-{{ sanitize(vtype) }}" role="tablist">
            {% for mtype in mtype_data.keys() %}
                <li class="nav-item">
                    <button class="nav-link {% if loop.first %}active{% endif %}"
                            id="tab-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}"
                            data-bs-toggle="tab"
                            data-bs-target="#pane-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}"
                            type="button">{{ mtype|title }}</button>
                </li>
            {% endfor %}
            </ul>

            <div class="tab-content" id="mtypeContent-{{ sanitize(vtype) }}">
            {% for mtype, groups in mtype_data.items() %}
                <div class="tab-pane fade {% if loop.first %}show active{% endif %}"
                     id="pane-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}" role="tabpanel">

                    <div id="mtype-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}"></div>

                    {% for group_name, subgroups in groups.items() %}
                    <div class="group-block"
                         id="group-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(group_name) }}">
                        <h3 class="group-toggle" onclick="toggleGroup(this)">{{ group_name }}</h3>
                        <div class="group-body" style="display:none;">

                        {% for sg_key, sg_data in subgroups.items() %}
                          {% if sg_key != "" %}
                          {# --- Subgroup accordeon (e.g. Chr 21) --- #}
                          <div class="subgroup-block"
                               id="sg-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(group_name) }}-{{ sanitize(sg_key) }}">
                              <h4 class="subgroup-toggle" onclick="toggleSubgroup(this)">{{ sg_data.name }}</h4>
                              <div class="subgroup-body" style="display:none;">
                              {# Sub-subgroups (e.g. "All" inside Chr 21) #}
                              {% for ssg_key, ssg_data in sg_data.subgroups.items() %}
                              <div class="subsubgroup-block"
                                   id="ssg-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sg_key) }}-{{ sanitize(ssg_key) }}">
                                  <h5 class="subsubgroup-toggle" onclick="toggleSubsubgroup(this)">{{ ssg_data.name }}</h5>
                                  <div class="subsubgroup-body" style="display:none;">
                                  {% for sec_key, sec_data in ssg_data.sections.items() %}
                                  <div class="section-card"
                                       id="sec-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sec_key) }}">
                                      <h4 class="section-toggle" onclick="toggleSection(this)">{{ sec_data.name }}</h4>
                                      <div class="section-body" style="display:none;">{{ sec_data.content|safe }}</div>
                                  </div>
                                  {% endfor %}
                                  </div>
                              </div>
                              {% endfor %}
                              {# Direct sections (ptype ones, not in sub-subgroup) #}
                              {% for sec_key, sec_data in sg_data.sections.items() %}
                              <div class="section-card"
                                   id="sec-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sec_key) }}">
                                  <h4 class="section-toggle" onclick="toggleSection(this)">{{ sec_data.name }}</h4>
                                  <div class="section-body" style="display:none;">{{ sec_data.content|safe }}</div>
                              </div>
                              {% endfor %}
                              </div>
                          </div>
                          {% else %}
                          {# --- Sections directly (no subgroup) --- #}
                          {% for sec_key, sec_data in sg_data.sections.items() %}
                          <div class="section-card"
                               id="sec-{{ sanitize(vtype) }}-{{ sanitize(mtype) }}-{{ sanitize(sec_key) }}">
                              <h4 class="section-toggle" onclick="toggleSection(this)">{{ sec_data.name }}</h4>
                              <div class="section-body" style="display:none;">{{ sec_data.content|safe }}</div>
                          </div>
                          {% endfor %}
                          {% endif %}
                        {% endfor %}

                        </div>
                    </div>
                    {% endfor %}

                </div>
            {% endfor %}
            </div>

        </div>
    {% endfor %}
    </div>

</div>

<!-- Image zoom modal -->
<div id="img-modal" onclick="this.style.display='none'">
    <img id="modal-img" src="" alt="Zoomed">
</div>

<script>
// Group toggle
function toggleGroup(el) {
    el.classList.toggle('open');
    var body = el.nextElementSibling;
    body.style.display = body.style.display === 'none' ? 'block' : 'none';
}
// Subgroup toggle
function toggleSubgroup(el) {
    el.classList.toggle('open');
    var body = el.nextElementSibling;
    body.style.display = body.style.display === 'none' ? 'block' : 'none';
}
// Sub-subgroup toggle
function toggleSubsubgroup(el) {
    el.classList.toggle('open');
    var body = el.nextElementSibling;
    body.style.display = body.style.display === 'none' ? 'block' : 'none';
}
// Section toggle
function toggleSection(el) {
    el.classList.toggle('open');
    var body = el.nextElementSibling;
    body.style.display = body.style.display === 'none' ? 'block' : 'none';
}

// Image zoom
document.addEventListener('click', function(e) {
    if (e.target.classList.contains('zoomable')) {
        document.getElementById('modal-img').src = e.target.src;
        document.getElementById('img-modal').style.display = 'flex';
        document.getElementById('img-modal').style.alignItems = 'center';
        document.getElementById('img-modal').style.justifyContent = 'center';
    }
});

// Smooth scroll for sidebar links
document.querySelectorAll('#sidebar .nav-link').forEach(function(link) {
    link.addEventListener('click', function(e) {
        var href = this.getAttribute('href');
        if (href && href.startsWith('#')) {
            e.preventDefault();
            var target = document.getElementById(href.substring(1));
            if (!target) return;

            // Extract vtype / mtype from the anchor prefix
            var secMatch  = href.match(/^#sec-([^-]+)-([^-]+)-/);
            var sgMatch   = href.match(/^#sg-([^-]+)-([^-]+)-/);
            var ssgMatch  = href.match(/^#ssg-([^-]+)-([^-]+)-/);
            var vtypeMatch = href.match(/^#vtype-([^-]+)/);
            var mtypeMatch = href.match(/^#mtype-([^-]+)-([^-]+)/);

            var vtypeId = null, mtypeId = null;
            if (secMatch)  { vtypeId = secMatch[1];  mtypeId = secMatch[2];  }
            else if (sgMatch)   { vtypeId = sgMatch[1];   mtypeId = sgMatch[2];   }
            else if (ssgMatch)  { vtypeId = ssgMatch[1];  mtypeId = ssgMatch[2];  }
            else if (vtypeMatch){ vtypeId = vtypeMatch[1]; }
            else if (mtypeMatch){ vtypeId = mtypeMatch[1]; mtypeId = mtypeMatch[2]; }

            // Activate correct Bootstrap tabs
            if (vtypeId) {
                var tabEl = document.getElementById('tab-' + vtypeId);
                if (tabEl) new bootstrap.Tab(tabEl).show();
            }
            if (vtypeId && mtypeId) {
                var tabEl2 = document.getElementById('tab-' + vtypeId + '-' + mtypeId);
                if (tabEl2) new bootstrap.Tab(tabEl2).show();
            }

            // Open every ancestor accordion body that is currently closed
            function openAncestors(el) {
                var node = el.parentElement;
                while (node && node.id !== 'main-content') {
                    if (node.classList.contains('group-body') ||
                        node.classList.contains('subgroup-body') ||
                        node.classList.contains('subsubgroup-body') ||
                        node.classList.contains('section-body')) {
                        if (node.style.display === 'none') {
                            node.style.display = 'block';
                            var toggle = node.previousElementSibling;
                            if (toggle) toggle.classList.add('open');
                        }
                    }
                    node = node.parentElement;
                }
            }
            openAncestors(target);

            // If the target is a section-card, also open its own content
            if (target.classList.contains('section-card')) {
                var secBody = target.querySelector(':scope > .section-body');
                if (secBody && secBody.style.display === 'none') {
                    secBody.style.display = 'block';
                    var secToggle = target.querySelector(':scope > .section-toggle');
                    if (secToggle) secToggle.classList.add('open');
                }
            }

            setTimeout(function() {
                target.scrollIntoView({behavior: 'smooth', block: 'start'});
            }, 300);
        }
    });
});

// Open first group in each mtype tab; close all others
document.querySelectorAll('.tab-pane').forEach(function(pane) {
    pane.querySelectorAll('.group-toggle').forEach(function(el, idx) {
        if (idx === 0) {
            el.classList.add('open');
            el.nextElementSibling.style.display = 'block';
            // Also open first subgroup (if any) inside first group
            var firstSg = el.nextElementSibling.querySelector('.subgroup-toggle');
            if (firstSg) {
                firstSg.classList.add('open');
                firstSg.nextElementSibling.style.display = 'block';
            }
            // Also open first section inside first group (or subgroup)
            var firstSec = el.nextElementSibling.querySelector('.section-toggle');
            if (firstSec) {
                firstSec.classList.add('open');
                firstSec.nextElementSibling.style.display = 'block';
            }
        }
    });
});

// Sidebar collapsible links with hidden children
(function() {
    var allLinks = Array.from(document.querySelectorAll('#sidebar .nav-link[data-mtype]'));
    
    // Mark links that have hidden children
    allLinks.forEach(function(link, idx) {
        var mtype = link.getAttribute('data-mtype');
        var level = null;
        for (var cls of link.classList) {
            if (cls.match(/^level-/)) {
                level = cls;
                break;
            }
        }
        if (!level) return;
        
        // Check if next sibling links are hidden children
        var hasHidden = false;
        for (var i = idx + 1; i < allLinks.length; i++) {
            var nextLink = allLinks[i];
            var nextLevel = null;
            for (var cls of nextLink.classList) {
                if (cls.match(/^level-/)) {
                    nextLevel = cls;
                    break;
                }
            }
            if (!nextLevel) continue;
            
            var levelNum = parseInt(level.split('-')[1]);
            var nextLevelNum = parseInt(nextLevel.split('-')[1]);
            
            // If we hit a same-or-lower level, stop
            if (nextLevelNum <= levelNum) break;
            
            // If child is collapsed-by-default, mark parent
            if (nextLink.classList.contains('sidebar-collapsed')) {
                hasHidden = true;
                break;
            }
        }
        
        if (hasHidden) {
            link.classList.add('has-hidden-children');
        }
    });
    
    // Toggle children on click
    allLinks.forEach(function(link, idx) {
        if (!link.classList.contains('has-hidden-children')) return;
        
        link.addEventListener('click', function(e) {
            var mtype = this.getAttribute('data-mtype');
            var level = null;
            for (var cls of this.classList) {
                if (cls.match(/^level-/)) {
                    level = cls;
                    break;
                }
            }
            if (!level) return;
            
            var levelNum = parseInt(level.split('-')[1]);
            var isExpanded = this.classList.contains('expanded');
            
            // Toggle children visibility
            for (var i = idx + 1; i < allLinks.length; i++) {
                var child = allLinks[i];
                var childLevel = null;
                for (var cls of child.classList) {
                    if (cls.match(/^level-/)) {
                        childLevel = cls;
                        break;
                    }
                }
                if (!childLevel) continue;
                
                var childLevelNum = parseInt(childLevel.split('-')[1]);
                if (childLevelNum <= levelNum) break;  // Stop at sibling/parent
                
                // Only toggle immediate children (levelNum + 1)
                if (childLevelNum === levelNum + 1) {
                    if (isExpanded) {
                        child.classList.add('sidebar-collapsed');
                        child.classList.remove('expanded');
                        // Hide all deeper descendants
                        for (var j = i + 1; j < allLinks.length; j++) {
                            var desc = allLinks[j];
                            var descLevel = null;
                            for (var cls of desc.classList) {
                                if (cls.match(/^level-/)) {
                                    descLevel = cls;
                                    break;
                                }
                            }
                            if (!descLevel) continue;
                            var descLevelNum = parseInt(descLevel.split('-')[1]);
                            if (descLevelNum <= levelNum) break;
                            if (descLevelNum > childLevelNum) {
                                desc.classList.add('sidebar-collapsed');
                            }
                        }
                    } else {
                        child.classList.remove('sidebar-collapsed');
                    }
                }
            }
            
            this.classList.toggle('expanded');
            e.preventDefault();
        });
    });
})();
</script>

</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate barometer Biomarker HTML Report")
    parser.add_argument("-r", "--results", default="barometer_results", help="Results directory from barometer_analyze.py")
    parser.add_argument("-o", "--output", default="barometer_report.html", help="Output HTML report file")
    parser.add_argument("-e", "--embed-images", action="store_true", default=False,
                        help="Embed images as base64 (larger file but self-contained)")
    args = parser.parse_args()
    # Pass embed flag as module-level so helpers can use it
    global EMBED_IMAGES, RESULTS_DIR
    EMBED_IMAGES = args.embed_images
    RESULTS_DIR = args.results

    if not os.path.isdir(args.results):
        log.error(f"Results directory not found: {args.results}")
        sys.exit(1)

    log.info(f"Building report from {args.results}...")
    report, manifest, global_content = build_report_data(args.results)

    log.info(f"Rendering HTML...")
    template = Template(HTML_TEMPLATE)
    html = template.render(
        report=report,
        manifest=manifest,
        global_content=global_content,
        sanitize=sanitize_id,
    )

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(html)

    log.info(f"Report written to {args.output}")
    log.info(f"Open in browser to view the interactive report.")


if __name__ == "__main__":
    main()
