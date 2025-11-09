#!/usr/bin/env python3
"""Test savgol_filter numerical consistency across scipy versions.

This script tests whether scipy.signal.savgol_filter produces identical
results between scipy 1.4.1 (Jan 2020) and scipy 1.16.2 (2025).

The test uses EXACT parameters from the virtual screen pipeline:
- window_length=5
- polyorder=3
- axis=1 (smooth each row independently)

Usage:
    # Modern scipy (1.16.2)
    pixi run python scripts/test_savgol_versions.py

    # Old scipy (1.4.1)
    pixi run --environment scipy14 python scripts/test_savgol_versions.py
"""

import numpy as np
import scipy
from scipy.signal import savgol_filter

# Fixed seed for reproducibility across runs
np.random.seed(42)

# Generate test data matching real pipeline dimensions
n_rows = 100
n_cols = 12  # 12 radial bins (features 5-16)
test_data = np.random.randn(n_rows, n_cols)

# Apply savgol_filter with EXACT pipeline parameters
# (from haghighi_mito/vectorized_slope.py:86)
smoothed = savgol_filter(test_data, window_length=5, polyorder=3, axis=1)

# Print environment info
print("=" * 80)
print("ENVIRONMENT")
print("=" * 80)
print(f"scipy version: {scipy.__version__}")
print(f"numpy version: {np.__version__}")

print()
print("=" * 80)
print("SAVGOL_FILTER TEST")
print("=" * 80)
print(f"Input shape: {test_data.shape}")
print(f"Input checksum: {np.sum(test_data):.15f}")
print(f"Output checksum: {np.sum(smoothed):.15f}")
print(f"Output std: {np.std(smoothed):.15f}")
print(f"Output mean: {np.mean(smoothed):.15f}")

print()
print("=" * 80)
print("SAMPLE VALUES (first 3 rows, first 6 columns)")
print("=" * 80)
for i in range(3):
    print(f"Row {i}: {smoothed[i][:6]}")

print()
print("=" * 80)
print("HASH TEST (for exact byte-level comparison)")
print("=" * 80)
# Hash first row for exact comparison
print(f"Row 0 hash: {hash(smoothed[0].tobytes())}")
print(f"Row 0 values: {smoothed[0]}")
