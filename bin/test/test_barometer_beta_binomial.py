import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import barometer_analyze


class BetaBinomialBridgeTest(unittest.TestCase):
    def test_technical_replicates_and_model_results_are_mapped(self):
        df = pd.DataFrame({
            "control::bio1::tech1::espf::successes": [1],
            "control::bio1::tech1::espf::trials": [10],
            "control::bio1::tech2::espf::successes": [2],
            "control::bio1::tech2::espf::trials": [20],
            "case::bio2::rep1::espf::successes": [8],
            "case::bio2::rep1::espf::trials": [40],
        }, index=[7])
        sample_info = [
            {"col": "control::bio1::tech1::espf", "group": "control", "sample": "bio1", "rep": "tech1"},
            {"col": "control::bio1::tech2::espf", "group": "control", "sample": "bio1", "rep": "tech2"},
            {"col": "case::bio2::rep1::espf", "group": "case", "sample": "bio2", "rep": "rep1"},
        ]

        def fake_rscript(command, **kwargs):
            counts = pd.read_csv(command[2], sep="\t")
            self.assertEqual(len(counts), 3)
            self.assertEqual(counts.loc[counts["group"] == "control", "successes"].sum(), 3)
            self.assertEqual(counts.loc[counts["group"] == "control", "trials"].sum(), 30)
            results = pd.DataFrame([
                {"row_id": "7", "test": "global", "group1": "", "group2": "", "statistic": 2.0,
                 "estimate": None, "mean1": None, "mean2": None, "p_value": 0.15, "status": "ok", "dispersion": 0.1},
                {"row_id": "7", "test": "mean", "group1": "case", "group2": "", "statistic": None,
                 "estimate": None, "mean1": 0.2, "mean2": None, "p_value": None, "status": "ok", "dispersion": 0.1},
                {"row_id": "7", "test": "mean", "group1": "control", "group2": "", "statistic": None,
                 "estimate": None, "mean1": 0.1, "mean2": None, "p_value": None, "status": "ok", "dispersion": 0.1},
                {"row_id": "7", "test": "pairwise", "group1": "case", "group2": "control", "statistic": None,
                 "estimate": 0.8, "mean1": 0.2, "mean2": 0.1, "p_value": 0.03, "status": "ok", "dispersion": 0.1},
            ])
            results.to_csv(command[3], sep="\t", index=False)
            return None

        with patch.object(barometer_analyze.subprocess, "run", side_effect=fake_rscript):
            result = barometer_analyze.run_beta_binomial_analysis(
                df, sample_info, ["case", "control"], tempfile.gettempdir()
            )

        self.assertEqual(result.loc[7, "beta_binomial_status"], "ok")
        self.assertEqual(result.loc[7, "beta_binomial_stat"], 2.0)
        self.assertEqual(result.loc[7, "beta_binomial_mean_case"], 0.2)
        self.assertEqual(result.loc[7, "beta_binomial_diff_case_vs_control"], 0.1)
        self.assertEqual(result.loc[7, "beta_binomial_pval_case_vs_control"], 0.03)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is available in the Barometer container")
    def test_glmmtmb_detects_synthetic_condition_effect(self):
        rows = []
        for group, counts in (("control", [1, 5, 8, 3, 12, 6]), ("case", [30, 45, 55, 25, 60, 40])):
            for index, successes in enumerate(counts, start=1):
                rows.append({
                    "row_id": "feature1",
                    "group": group,
                    "sample": f"{group}_{index}",
                    "replicate": "rep1",
                    "successes": successes,
                    "trials": 100,
                })
                rows.append({
                    "row_id": "allzero",
                    "group": group,
                    "sample": f"{group}_{index}",
                    "replicate": "rep1",
                    "successes": 0,
                    "trials": 100,
                })

        script_path = Path(barometer_analyze.__file__).with_name("barometer_beta_binomial.R")
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "counts.tsv"
            output_path = Path(temp_dir) / "results.tsv"
            pd.DataFrame(rows).to_csv(input_path, sep="\t", index=False)
            subprocess.run(
                ["Rscript", str(script_path), str(input_path), str(output_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            results = pd.read_csv(output_path, sep="\t")

        global_result = results.loc[results["test"] == "global"].iloc[0]
        pairwise_results = results.loc[results["test"] == "pairwise"]
        self.assertFalse(pairwise_results.empty, results.to_string(index=False))
        pairwise_result = pairwise_results.iloc[0]
        self.assertEqual(global_result["status"], "ok")
        self.assertLess(global_result["p_value"], 0.05)
        self.assertEqual(pairwise_result["group1"], "case")
        self.assertEqual(pairwise_result["group2"], "control")
        self.assertGreater(pairwise_result["estimate"], 0)
        self.assertGreater(pairwise_result["mean1"], pairwise_result["mean2"])
        zero_result = results.loc[
            (results["row_id"] == "allzero") & (results["test"] == "global")
        ].iloc[0]
        self.assertEqual(zero_result["status"], "no_outcome_variation")


if __name__ == "__main__":
    unittest.main()