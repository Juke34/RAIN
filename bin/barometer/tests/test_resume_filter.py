#!/usr/bin/env python3
"""Test the --resume mode: TSV loading + task filtering.

This is a *synthetic* test. It does not import the ``barometer`` package
(which would pull in seaborn / sklearn / statsmodels). Instead it replicates
the two pieces of resume logic that live in ``barometer.__main__.main()``:

  1. Loading ``failed_analyses.tsv`` into a set of
     ``(vtype, mtype, section_key)`` string tuples.
  2. Filtering the per-vtype ``tasks`` list so that only the sections listed
     in the TSV are re-run.

It asserts that, given a TSV listing one failed section, exactly that section
is kept and the others are dropped.

Run:
    python3 tests/test_resume_filter.py
"""

import os
import tempfile

import pandas as pd


def load_resume_keys(path):
    """Replica of the --resume TSV loading in __main__.main()."""
    if not os.path.isfile(path):
        raise SystemExit(f"--resume file not found: {path}")
    df_failed = pd.read_csv(path, sep="\t", low_memory=False)
    required_cols = {"vtype", "mtype", "section_key"}
    missing = required_cols - set(df_failed.columns)
    if missing:
        raise SystemExit(f"--resume file {path} is missing columns: {sorted(missing)}")
    return set(zip(df_failed["vtype"].astype(str),
                   df_failed["mtype"].astype(str),
                   df_failed["section_key"].astype(str)))


def filter_tasks(vtype, tasks, resume_keys):
    """Replica of the per-vtype task filter in __main__.main().

    Each task is a 9-tuple whose element 5 is the key ``(mtype, section_key)``.
    """
    return [t for t in tasks if (vtype, str(t[5][0]), str(t[5][1])) in resume_keys]


def main():
    # Build a failed_analyses.tsv with one failed section (Task B).
    with tempfile.TemporaryDirectory() as tmp:
        tsv = os.path.join(tmp, "failed_analyses.tsv")
        pd.DataFrame([{
            "vtype": "espf", "mtype": "feature", "section_key": "secB",
            "section_name": "Task B", "section_outdir": "/x/secB",
            "reason": "timeout", "detail": "exceeded 2400s",
            "timestamp": "2026-01-01 00:00:00",
        }]).to_csv(tsv, sep="\t", index=False)

        resume_keys = load_resume_keys(tsv)
        assert resume_keys == {("espf", "feature", "secB")}, resume_keys

        # Simulate the tasks built for vtype='espf' (9-tuples, key at index 5).
        vtype = "espf"
        tasks = [
            ("df", ["c"], {}, "o", "Task A", ("feature", "secA"), "st", None, 10),
            ("df", ["c"], {}, "o", "Task B", ("feature", "secB"), "st", None, 10),
            ("df", ["c"], {}, "o", "Task C", ("feature", "secC"), "st", None, 10),
        ]
        filtered = filter_tasks(vtype, tasks, resume_keys)

        assert len(filtered) == 1, f"Expected 1 kept task, got {len(filtered)}"
        assert filtered[0][4] == "Task B", f"Expected Task B kept, got {filtered[0][4]}"

        # A vtype not present in the TSV keeps nothing.
        assert filter_tasks("espr", tasks, resume_keys) == []

    print("TEST PASSED: --resume keeps only the failed sections")


if __name__ == "__main__":
    main()
