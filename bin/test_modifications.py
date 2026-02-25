#!/usr/bin/env python3
"""Test script to validate the new by-parent-type aggregation modifications"""

import sys
sys.path.insert(0, '/Users/jacquesdainat/git/Juke34/rain/bin')

try:
    # Test imports
    from pluviometer.rain_file_writers import AggregateFileWriter
    from collections import defaultdict
    
    print("✓ Imports successful")
    
    # Test that the new method exists
    if hasattr(AggregateFileWriter, 'write_rows_with_data_by_parent_type'):
        print("✓ New method 'write_rows_with_data_by_parent_type' exists in AggregateFileWriter")
    else:
        print("✗ Method 'write_rows_with_data_by_parent_type' NOT FOUND")
        sys.exit(1)
    
    # Test the method signature
    import inspect
    sig = inspect.signature(AggregateFileWriter.write_rows_with_data_by_parent_type)
    params = list(sig.parameters.keys())
    expected_params = ['self', 'record_id', 'parent_list', 'aggregate_id', 'aggregation_mode', 'counter_dict']
    
    if params == expected_params:
        print(f"✓ Method signature correct: {params}")
    else:
        print(f"✗ Method signature mismatch. Expected: {expected_params}, Got: {params}")
        sys.exit(1)
    
    print("\n✓ All validation tests passed!")
    print("\nThe modifications add:")
    print("  - Aggregation by (ParentType, AggregateType) at SeqID level")
    print("  - Aggregation by (ParentType, AggregateType) at genome level")
    print("  - New counters: *_by_parent_type for each aggregation mode")
    
except Exception as e:
    print(f"✗ Error during validation: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
