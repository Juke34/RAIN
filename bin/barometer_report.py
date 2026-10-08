#!/usr/bin/env python3
"""Compatibility shim for the legacy ``barometer_report.py`` entry point.

The implementation now lives in the ``barometer`` package (see
``bin/barometer/report.py``). This module is kept so that existing
invocations of ``python barometer_report.py ...`` keep working. It simply
delegates to the package's ``main``.
"""

import os
import sys

# Make sure the package is importable when this file is run directly
# (e.g. ``python bin/barometer_report.py``) without bin/ on sys.path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from barometer.report import main  # noqa: E402

if __name__ == "__main__":
    main()

