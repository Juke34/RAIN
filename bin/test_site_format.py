"""
Tests for the per-site pluviometer output format (standard 16-column format)
and its downstream processing by drip.py.

Run with: python -m pytest bin/test_site_format.py -v
"""

import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

# pluviometer.utils imports Bio.SeqFeature at module level (only for type
# annotations on location helpers, unused by RNASiteVariantData). Stub it so
# the unit tests run without BioPython installed.
if "Bio" not in sys.modules:
    bio = types.ModuleType("Bio")
    seqfeature = types.ModuleType("Bio.SeqFeature")
    seqfeature.SimpleLocation = object
    seqfeature.CompoundLocation = object
    bio.SeqFeature = seqfeature
    sys.modules["Bio"] = bio
    sys.modules["Bio.SeqFeature"] = seqfeature

from pluviometer.site_file_writer import SiteFileWriter, SITE_HEADER
from pluviometer.utils import RNASiteVariantData

HEADER = "\t".join(SITE_HEADER)


def make_site(seqid="chr1", position=9, reference=0, strand=1,
              counts=(60, 0, 40, 0), coverage=100) -> RNASiteVariantData:
    return RNASiteVariantData(
        seqid=seqid,
        position=position,
        reference=reference,
        strand=strand,
        coverage=coverage,
        mean_quality=0.0,
        frequencies=np.array([*counts, 0], dtype=np.int64),
        score=0.0,
    )


def test_site_writer_standard_format(tmp_path: Path) -> None:
    out = tmp_path / "sites.tsv"
    with out.open("w") as handle:
        writer = SiteFileWriter(handle, cov_threshold=30, edit_threshold=3)
        writer.write_header()
        # Site in a feature: ref A, 60 A / 40 G reads
        writer.write_site(make_site(), ["geneA"])
        # Site outside any feature: ref G, 50 G / 10 A reads
        writer.write_site(make_site(position=29, reference=2, counts=(10, 0, 50, 0)), [])
        # Below coverage threshold: skipped
        writer.write_site(make_site(position=39, coverage=10), [])

    lines = out.read_text().splitlines()
    assert lines[0] == HEADER
    assert len(lines) == 3  # header + 2 sites

    row = lines[1].split("\t")
    assert len(row) == 16
    seqid, parent_ids, sid, mtype, ptype, ftype, ctype, mode, start, end, strand, total, obs, qual, sbp, rbp = row
    assert seqid == "chr1"
    assert parent_ids == ".,geneA"
    assert sid == "site:chr1:10:1"
    assert mtype == "site" and ptype == "." and ftype == "site" and ctype == "."
    assert mode == "all_sites"
    assert start == "10" and end == "10" and strand == "1"
    assert total == "1"
    # One-hot on the reference base A
    assert obs == "1,0,0,0" and qual == "1,0,0,0"
    # Site pairings: single 1 on AA (index 0)
    assert sbp == "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
    # Read pairings: row A = (60, 0, 40, 0)
    assert rbp == "60,0,40,0,0,0,0,0,0,0,0,0,0,0,0,0"

    row2 = lines[2].split("\t")
    assert row2[1] == "."  # outside any feature
    assert row2[2] == "site:chr1:30:1"
    assert row2[12] == "0,0,1,0"  # one-hot on G
    assert row2[14] == "0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0"  # GG at index 10 (2*4+2)
    assert row2[15] == "0,0,0,0,0,0,0,0,10,0,50,0,0,0,0,0"  # row G = (10, 0, 50, 0)


def test_site_writer_edit_threshold(tmp_path: Path) -> None:
    out = tmp_path / "sites.tsv"
    with out.open("w") as handle:
        writer = SiteFileWriter(handle, cov_threshold=30, edit_threshold=3)
        writer.write_header()
        # Non-ref count 2 < edit_threshold 3 → zeroed
        writer.write_site(make_site(counts=(98, 0, 2, 0)), [])

    row = out.read_text().splitlines()[1].split("\t")
    assert row[15] == "98,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"


def test_drip_on_site_files(tmp_path: Path) -> None:
    # Site 10: A-to-G edited in 2 samples; site 30: A-to-G edited in 1 sample only
    def write_sites(path: Path, rows: list[str]) -> None:
        path.write_text(HEADER + "\n" + "".join(r + "\n" for r in rows))

    s1 = tmp_path / "s1.tsv"
    s2 = tmp_path / "s2.tsv"
    s3 = tmp_path / "s3.tsv"
    write_sites(s1, [
        "chr1\t.,geneA\tsite:chr1:10:1\tsite\t.\tsite\t.\tall_sites\t10\t10\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t60,0,40,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "chr1\t.\tsite:chr1:30:1\tsite\t.\tsite\t.\tall_sites\t30\t30\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t50,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0",
    ])
    write_sites(s2, [
        "chr1\t.,geneA\tsite:chr1:10:1\tsite\t.\tsite\t.\tall_sites\t10\t10\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t70,0,30,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "chr1\t.\tsite:chr1:30:1\tsite\t.\tsite\t.\tall_sites\t30\t30\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t40,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0",
    ])
    write_sites(s3, [
        "chr1\t.,geneA\tsite:chr1:10:1\tsite\t.\tsite\t.\tall_sites\t10\t10\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "chr1\t.\tsite:chr1:30:1\tsite\t.\tsite\t.\tall_sites\t30\t30\t1\t1\t1,0,0,0\t1,0,0,0\t1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\t60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
    ])

    import os
    script = Path(__file__).parent / "drip.py"
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        subprocess.run(
            [sys.executable, str(script), "--output", "drip_sites",
             f"{s1}:ctl:s1:rep1", f"{s2}:ctl:s2:rep1", f"{s3}:trt:s3:rep1"],
            check=True,
        )
    finally:
        os.chdir(cwd)

    df = pd.read_csv(tmp_path / "drip_sites_espr" / "drip_sites_espr_AG.tsv", sep="\t")
    assert sorted(df["Start"]) == [10, 30]
    row = df[df["Start"] == 10].iloc[0]
    assert row["ctl::s1::rep1::espr::successes"] == 40
    assert row["ctl::s1::rep1::espr::trials"] == 100
    assert row["ctl::s2::rep1::espr::successes"] == 30
    assert row["trt::s3::rep1::espr::successes"] == 0
    assert row["trt::s3::rep1::espr::trials"] == 100
    assert row["ParentIDs"] == ".,geneA"
    assert df[df["Start"] == 30].iloc[0]["ParentIDs"] == "."
