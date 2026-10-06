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
            report_non_qualified=True,
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
        output_prefix = str(Path(self.temp_dir.name) / "lowcov")
        drip.merge_samples(
            {str(self.input_path): ("control", "sample1", "rep1")},
            output_prefix,
            min_cov=11,
            report_non_qualified=True,
        )

        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertTrue(pd.isna(espf.loc[0, "control::sample1::rep1::espf::successes"]))
        self.assertTrue(pd.isna(espf.loc[0, "control::sample1::rep1::espf::trials"]))

    def test_preserve_covered_rows_keeps_zero_successes(self):
        row = pd.read_csv(self.input_path, sep="\t").iloc[0].to_dict()
        row["SiteBasePairingsQualified"] = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
        zero_input = Path(self.temp_dir.name) / "zero.tsv"
        pd.DataFrame([row]).to_csv(zero_input, sep="\t", index=False)
        output_prefix = str(Path(self.temp_dir.name) / "covered")

        drip.merge_samples(
            {str(zero_input): ("control", "sample1", "rep1")},
            output_prefix,
            min_samples_pct=100,
            min_group_pct=100,
            preserve_covered_features=True,
        )

        espf = pd.read_csv(f"{output_prefix}_espf_AC.tsv", sep="\t")
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf"], 0.0)
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::successes"], 0)
        self.assertEqual(espf.loc[0, "control::sample1::rep1::espf::trials"], 10)


if __name__ == "__main__":
    unittest.main()