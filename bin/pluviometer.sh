#!/bin/bash
# Wrapper script to call pluviometer module from bin directory
# This maintains backward compatibility when pluviometer.py was in bin/

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run pluviometer as a Python module
exec python -m pluviometer "$@"
