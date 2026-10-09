import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import drip


class DripCountOutputTest(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.input_path = Path(self.temp_dir.name) / "sample.tsv"
        row = {
            "SeqID": "chr1",
            "ParentIDs": ".",
            "ID": "gene1",
            "Mtype": "aggregate",
            "Ptype": ".",
            "Type": "feature",
            "Ctype": ".",
            "Mode": ".",
            "Start": "1",
            "End": "10",
            "Strand": "+",
            "TotalSites": "10",
            "ObservedBases": "10,2,0,0",
            "QualifiedBases": "10,2,0,0",
            "SiteBasePairingsQualified": "0,3,2,0,0,0,0,1,0,0,0,0,0,0,0,0",
            "ReadBasePairingsQualified": "10,3,2,1,4,2,3,1,0,0,0,0,0,0,0,0",
        }
        pd.DataFrame([row]).to_csv(self.input_path, sep="\t", index=False)

    def test_merge_keeps_ratios_and_integer_counts_in_same_table(self):
        output_prefix = str(Path(self.temp_dir.name) / "drip")
        drip.merge_samples(
            {str(self.input_path): ("control", "sample1", "rep1")},
            output_prefix,
            preserve_covered_features=True,
        )

        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        espr = pd.read_csv(f"{output_prefix}_espr_AC.tsv", sep="\t")

        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf"], 0.3)
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::successes"], 3)
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::trials"], 10)
        self.assertEqual(espr.loc[0, "control::sample1::rep1::espr"], 0.1875)
        self.assertEqual(espr.loc[0, "control::sample1::rep1::espr::successes"], 3)
        self.assertEqual(espr.loc[0, "control::sample1::rep1::espr::trials"], 16)

        ag = pd.read_csv(f"{output_prefix}_espf_AG.tsv", sep="\t")
        ct = pd.read_csv(f"{output_prefix}_espr_CT.tsv", sep="\t")
        self.assertEqual(ag.loc[0, "control::sample1::rep1::espf::successes"], 2)
        self.assertEqual(ag.loc[0, "control::sample1::rep1::espf::trials"], 10)
        self.assertEqual(ct.loc[0, "control::sample1::rep1::espr::successes"], 1)
        self.assertEqual(ct.loc[0, "control::sample1::rep1::espr::trials"], 10)

    def test_low_coverage_counts_are_missing(self):
        # min_cov=11 puts every cell below the coverage threshold → all NA.
        # With no coverage filter flag set, the row is kept (all rows pass)
        # but its count columns are NA.
        output_prefix = str(Path(self.temp_dir.name) / "lowcov")
        drip.merge_samples(
            {str(self.input_path): ("control", "sample1", "rep1")},
            output_prefix,
            min_cov=11,
            preserve_covered_features=True,
        )

        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)
        self.assertTrue(pd.isna(espf.loc[0, "control::sample1::rep1::espf::successes"]))
        self.assertTrue(pd.isna(espf.loc[0, "control::sample1::rep1::espf::trials"]))

    def test_preserve_covered_features_only_adds_count_columns(self):
        # With the flag OFF, count columns are not written to the output.
        output_prefix = str(Path(self.temp_dir.name) / "no_counts")
        drip.merge_samples(
            {str(self.input_path): ("control", "sample1", "rep1")},
            output_prefix,
        )
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertNotIn("control::sample1::rep1::espf::successes", espf.columns)
        self.assertNotIn("control::sample1::rep1::espf::trials", espf.columns)

        # With the flag ON, count columns are written.
        output_prefix = str(Path(self.temp_dir.name) / "with_counts")
        drip.merge_samples(
            {str(self.input_path): ("control", "sample1", "rep1")},
            output_prefix,
            preserve_covered_features=True,
        )
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::successes"], 3)
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::trials"], 10)

    def test_preserve_covered_features_does_not_bypass_row_filter(self):
        # A row with zero editing in every sample (covered but no edit) must be
        # filtered out even when --preserve-covered-features is set: the flag
        # only adds count columns, it does not relax row filtering.
        row = pd.read_csv(self.input_path, sep="\t").iloc[0].to_dict()
        row["SiteBasePairingsQualified"] = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        zero_input = Path(self.temp_dir.name) / "zero.tsv"
        pd.DataFrame([row]).to_csv(zero_input, sep="\t", index=False)
        output_prefix = str(Path(self.temp_dir.name) / "zero_filtered")

        drip.merge_samples(
            {str(zero_input): ("control", "sample1", "rep1")},
            output_prefix,
            min_samples_pct=100,
            min_group_samples_pct=100,
            min_group_samples_edited=1,
            preserve_covered_features=True,
        )

        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            # No rows at all → file not written.
            pass

    def _make_no_edit_input(self):
        """Input file identical to self.input_path but with zero editing."""
        row = pd.read_csv(self.input_path, sep="\t").iloc[0].to_dict()
        row["SiteBasePairingsQualified"] = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        path = Path(self.temp_dir.name) / "no_edit.tsv"
        pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
        return path

    def _two_group_samples(self, no_edit_path):
        """4 samples in 2 groups: small (1, edited) / big (3, no edit).

        merge_samples keys its dict by file path, so each sample needs its
        own copy of the no-edit input.
        """
        big_paths = {}
        for name in ("b1", "b2", "b3"):
            p = Path(self.temp_dir.name) / f"no_edit_{name}.tsv"
            pd.read_csv(no_edit_path, sep="\t").to_csv(p, sep="\t", index=False)
            big_paths[name] = str(p)
        return {
            str(self.input_path): ("small", "s1", "rep1"),
            big_paths["b1"]: ("big", "b1", "rep1"),
            big_paths["b2"]: ("big", "b2", "rep1"),
            big_paths["b3"]: ("big", "b3", "rep1"),
        }

    def test_min_group_samples_and_pct_same_group(self):
        # Two groups: "small" (1 sample, edited) and "big" (3 samples, no edit
        # but covered).  Coverage filters count non-NA cells, so both groups
        # are 100% covered; the editing stage is disabled by default.
        no_edit_path = self._make_no_edit_input()
        samples = self._two_group_samples(no_edit_path)

        # min_group_samples_pct=75 alone: small group 1/1 = 100% >= 75% → row
        # kept (the permissive single-sample case).
        output_prefix = str(Path(self.temp_dir.name) / "pct_only")
        drip.merge_samples(samples, output_prefix, min_group_samples_pct=75)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        # min_group_samples_pct=75 + min_group_samples=2: small group has
        # only 1 covered sample < 2 → not qualified; big group is 100%
        # covered and 100% >= 75% → qualified.  Editing stage: small group
        # has 1 edited sample → row kept.
        output_prefix = str(Path(self.temp_dir.name) / "pct_and_n2")
        drip.merge_samples(
            samples, output_prefix, min_group_samples_pct=75, min_group_samples=2
        )
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        # min_group_samples_pct=75 + min_group_samples=4: no group has >= 4
        # covered samples (small: 1, big: 3) → coverage stage fails → row
        # filtered, even though the small group has an edited sample.
        output_prefix = str(Path(self.temp_dir.name) / "pct_and_n4")
        drip.merge_samples(
            samples, output_prefix, min_group_samples_pct=75, min_group_samples=4
        )
        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            pass  # No rows → file not written.

    def test_min_samples_and_pct_same_counter(self):
        # 2 samples total, both covered (1 edited, 1 covered-without-edit).
        # Coverage is 2/2 = 100%, so min_samples_pct=50 alone keeps the row.
        # The editing stage is disabled (min_group_samples_edited=0) to
        # isolate the coverage AND logic.
        no_edit_path = self._make_no_edit_input()
        samples = {
            str(self.input_path): ("g1", "s1", "rep1"),
            str(no_edit_path): ("g1", "s2", "rep1"),
        }

        output_prefix = str(Path(self.temp_dir.name) / "pct50_only")
        drip.merge_samples(
            samples, output_prefix, min_samples_pct=50, min_group_samples_edited=0
        )
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        # min_samples_pct=100 + min_samples=3: coverage 2/2 passes the pct
        # but 2 < 3 → filtered (AND on the same global counter).
        output_prefix = str(Path(self.temp_dir.name) / "pct100_n3")
        drip.merge_samples(
            samples, output_prefix,
            min_samples_pct=100, min_samples=3, min_group_samples_edited=0,
        )
        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            pass  # No rows → file not written.

        # min_samples=1 alone: 1 covered sample → row kept.
        output_prefix = str(Path(self.temp_dir.name) / "n1_only")
        drip.merge_samples(samples, output_prefix, min_samples=1)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

    def test_min_group_samples_alone(self):
        # Same setup: with min_group_samples=1 alone, the small group
        # qualifies (1 covered sample) → row kept.
        no_edit_path = self._make_no_edit_input()
        samples = self._two_group_samples(no_edit_path)
        output_prefix = str(Path(self.temp_dir.name) / "n1_only")
        drip.merge_samples(samples, output_prefix, min_group_samples=1)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

    def test_coverage_counts_zeros(self):
        # A row covered in every sample but with zero editing everywhere:
        # coverage stage passes (non-NA counts).  With the default editing
        # stage (disabled) the row is kept — "no editing" is reported.
        # With min_group_samples_edited=1 the row is filtered.
        no_edit_path = self._make_no_edit_input()
        samples = {str(no_edit_path): ("g1", "s1", "rep1")}

        output_prefix = str(Path(self.temp_dir.name) / "cov_default_edit")
        drip.merge_samples(samples, output_prefix, min_samples_pct=100)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)
        self.assertEqual(espf.loc[0, "g1::s1::rep1::espf"], 0.0)

        output_prefix = str(Path(self.temp_dir.name) / "cov_edit1")
        drip.merge_samples(
            samples, output_prefix,
            min_samples_pct=100, min_group_samples_edited=1,
        )
        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            pass  # No rows → file not written.

    def test_min_group_samples_edited(self):
        # Two groups: "small" (1 sample, edited) and "big" (3 samples, no
        # edit but covered).  With the default (editing stage disabled) the
        # row is kept; min_group_samples_edited=1 also keeps it via the
        # small group.
        no_edit_path = self._make_no_edit_input()
        samples = self._two_group_samples(no_edit_path)

        output_prefix = str(Path(self.temp_dir.name) / "edited_default")
        drip.merge_samples(samples, output_prefix)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        output_prefix = str(Path(self.temp_dir.name) / "edited_n1")
        drip.merge_samples(samples, output_prefix, min_group_samples_edited=1)
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        # min_group_samples_edited=2: no group has >= 2 edited samples
        # (small has 1, big has 0) → row filtered.
        output_prefix = str(Path(self.temp_dir.name) / "edited_n2")
        drip.merge_samples(samples, output_prefix, min_group_samples_edited=2)
        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            pass  # No rows → file not written.

    def test_min_group_samples_pct_edited(self):
        # Two groups: "small" (1 sample, edited) and "big" (3 samples, no
        # edit but covered).  min_group_samples_pct_edited=50: small group
        # has 1/1 = 100% edited >= 50% → row kept.
        no_edit_path = self._make_no_edit_input()
        samples = self._two_group_samples(no_edit_path)

        output_prefix = str(Path(self.temp_dir.name) / "pct_edited_50")
        drip.merge_samples(
            samples, output_prefix, min_group_samples_pct_edited=50
        )
        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(len(espf), 1)

        # min_group_samples_pct_edited=50 + min_group_samples_edited=2:
        # small group has 1 edited sample < 2 → AND fails; big group has 0%
        # edited → fails. Row filtered.
        output_prefix = str(Path(self.temp_dir.name) / "pct_edited_50_n2")
        drip.merge_samples(
            samples, output_prefix,
            min_group_samples_pct_edited=50, min_group_samples_edited=2,
        )
        espf_path = f"{output_prefix}_espf_AC.tsv"
        if os.path.exists(espf_path):
            espf = pd.read_csv(espf_path, sep="\t")
            self.assertEqual(len(espf), 0)
        else:
            pass  # No rows → file not written.


if __name__ == "__main__":
    unittest.main()