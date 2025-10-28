#!/usr/bin/env python3
"""Test script to check GDAL WarpOptions parameters"""

from osgeo import gdal
import traceback

print(f"GDAL Version: {gdal.VersionInfo()}")

# Test 1: Check what parameters WarpOptions accepts
print("\nTesting GDAL WarpOptions parameters...")

try:
    # Test with cutline parameter
    options = gdal.WarpOptions(
        format='MEM',
        cutline='POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))'
    )
    print("✓ 'cutline' parameter works")
except Exception as e:
    print(f"✗ 'cutline' parameter failed: {e}")

try:
    # Test with cutlineWKT parameter (might be old version)
    options = gdal.WarpOptions(
        format='MEM',
        cutlineWKT='POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))'
    )
    print("✓ 'cutlineWKT' parameter works")
except Exception as e:
    print(f"✗ 'cutlineWKT' parameter failed: {e}")

try:
    # Test with cropToCutline parameter
    options = gdal.WarpOptions(
        format='MEM',
        cropToCutline=True
    )
    print("✓ 'cropToCutline' parameter works")
except Exception as e:
    print(f"✗ 'cropToCutline' parameter failed: {e}")

# Test 2: Check help documentation
print("\nChecking GDAL WarpOptions help:")
try:
    help_info = gdal.WarpOptions.__doc__
    if help_info:
        print(help_info[:500] + "..." if len(help_info) > 500 else help_info)
    else:
        print("No help documentation available")
except:
    print("Could not access help documentation")

# Test 3: Try creating a basic WarpOptions
print("\nTesting basic WarpOptions creation:")
try:
    basic_options = gdal.WarpOptions(format='MEM')
    print("✓ Basic WarpOptions creation works")
except Exception as e:
    print(f"✗ Basic WarpOptions creation failed: {e}")
    print(traceback.format_exc())