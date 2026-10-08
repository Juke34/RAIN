#!/usr/bin/env python
"""Entry point for the ``barometer`` package (analyze / merge / report).

The implementation lives in the ``barometer`` package (see
``bin/barometer/``). This wrapper dispatches to the right sub-command:

* ``--report ...``  → ``barometer.report.main()``   (HTML report)
* anything else     → ``barometer.__main__.main()`` (analyze / --merge)

It is also kept as a compatibility shim for the legacy
``barometer_analyze.py`` entry point.
"""

import os
import sys

# Limit implicit BLAS/LAPACK multi-threading so each of the -j worker
# processes uses at most 1 core (total CPU usage = n_jobs, not n_jobs^2).
# Must be set before numpy/scipy/sklearn are imported.
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')

# Make sure the package is importable when this file is run directly
# (e.g. ``python bin/barometer_wrapper.py``) without bin/ on sys.path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    if "--report" in sys.argv:
        # --report is a dispatch flag for this wrapper, not a real option of
        # the report parser → strip it before delegating.
        sys.argv.remove("--report")
        from barometer.report import main
    else:
        from barometer.__main__ import main
    main()
