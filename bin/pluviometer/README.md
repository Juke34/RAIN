# Pluviometer Package

## Description

Pluviometer is an RNA editing quantification tool that counts and analyzes editing sites from variant calling data in the context of genomic features.

## Structure

The pluviometer package contains the following modules:

- `__main__.py` - Main entry point of the program (formerly pluviometer.py)
- `rain_file_writers.py` - Output file writers (features and aggregates)
- `multi_counter.py` - Multi-level counters for editing sites
- `site_filter.py` - Site filtering based on thresholds
- `rna_site_variant_readers.py` - Readers for different formats (Reditools2/3, Jacusa2)
- `SeqFeature_extensions.py` - Custom extensions and methods for Bio.SeqFeature
- `utils.py` - Utility functions and data types

## Usage

From the `bin/` directory, you can call pluviometer in several ways:

### 1. As a Python module (recommended)
```bash
python -m pluviometer --sites SITES --gff GFF [OPTIONS]
```

### 2. Via the shell wrapper
```bash
./pluviometer.sh --sites SITES --gff GFF [OPTIONS]
```

### 3. Via the Python wrapper
```bash
python pluviometer_wrapper.py --sites SITES --gff GFF [OPTIONS]
```

All these methods are equivalent and maintain compatibility with the old `python pluviometer.py` call.

## Migration

The following files have been moved from `bin/` to `bin/pluviometer/`:
- pluviometer.py → pluviometer/__main__.py
- rain_file_writers.py → pluviometer/rain_file_writers.py
- multi_counter.py → pluviometer/multi_counter.py
- site_filter.py → pluviometer/site_filter.py
- rna_site_variant_readers.py → pluviometer/rna_site_variant_readers.py
- SeqFeature_extensions.py → pluviometer/SeqFeature_extensions.py
- utils.py → pluviometer/utils.py

All internal imports have been converted to relative imports to function as a Python package.
