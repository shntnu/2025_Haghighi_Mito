"""Tests for vectorized slope calculation functions."""

import numpy as np
import pytest

from haghighi_mito.vectorized_slope import find_end_slope2, find_end_slope2_vectorized


class TestVectorizedSlope:
    """Test suite for vectorized slope calculation."""

    def test_single_row_equivalence(self):
        """Test that vectorized version matches original for a single row."""
        # Create a simple profile with a clear peak
        data = np.array([0.1, 0.3, 0.5, 0.8, 1.0, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.08, 0.06, 0.05])

        # Original function
        last_peak_orig, slope_orig = find_end_slope2(data)

        # Vectorized function (reshape to 2D)
        result_vec = find_end_slope2_vectorized(data.reshape(1, -1))
        last_peak_vec, slope_vec = result_vec[0, 0], result_vec[0, 1]

        # Should match exactly
        assert last_peak_orig == last_peak_vec, f"Peak indices don't match: {last_peak_orig} vs {last_peak_vec}"
        assert np.isclose(slope_orig, slope_vec), f"Slopes don't match: {slope_orig} vs {slope_vec}"

    def test_multiple_rows_equivalence(self):
        """Test that vectorized version matches original for multiple rows."""
        # Generate random profiles
        np.random.seed(42)
        n_rows = 100
        n_cols = 16

        # Create profiles that look like radial distributions
        # (typically start low, peak in middle, end low)
        x = np.linspace(0, 2 * np.pi, n_cols)
        data_matrix = np.zeros((n_rows, n_cols))
        for i in range(n_rows):
            # Random sine-like curves with noise
            data_matrix[i] = np.sin(x + np.random.rand()) * np.random.rand() + np.random.randn(n_cols) * 0.1

        # Process with original function (row by row)
        results_orig = np.apply_along_axis(find_end_slope2, 1, data_matrix)

        # Process with vectorized function
        results_vec = find_end_slope2_vectorized(data_matrix)

        # Compare results
        np.testing.assert_array_equal(results_orig[:, 0], results_vec[:, 0], err_msg="Peak indices don't match")
        np.testing.assert_allclose(results_orig[:, 1], results_vec[:, 1], rtol=1e-10, err_msg="Slopes don't match")

    def test_edge_case_flat_profile(self):
        """Test handling of flat profiles with no peaks."""
        # Flat profile
        data = np.ones((5, 16))

        results = find_end_slope2_vectorized(data)

        # With a flat profile, the peak will be at the start or end
        # The function should handle this gracefully
        assert results.shape == (5, 2), "Output shape is incorrect"
        assert not np.any(np.isnan(results)), "Results contain NaN values"

    def test_edge_case_monotonic_increasing(self):
        """Test handling of monotonically increasing profiles."""
        # Monotonically increasing
        data = np.linspace(0, 1, 16).reshape(1, -1)

        result_orig = np.array(find_end_slope2(data[0]))
        result_vec = find_end_slope2_vectorized(data)

        np.testing.assert_array_equal(result_orig, result_vec[0])

    def test_edge_case_monotonic_decreasing(self):
        """Test handling of monotonically decreasing profiles."""
        # Monotonically decreasing
        data = np.linspace(1, 0, 16).reshape(1, -1)

        result_orig = np.array(find_end_slope2(data[0]))
        result_vec = find_end_slope2_vectorized(data)

        np.testing.assert_array_equal(result_orig, result_vec[0])

    def test_realistic_radial_distribution(self):
        """Test with profiles that resemble actual radial distributions."""
        # Typical radial distribution: starts low, peaks in middle, decreases
        n_bins = 16
        profiles = []

        # Control profile - relatively flat
        profiles.append(np.ones(n_bins) * 0.5 + np.random.randn(n_bins) * 0.05)

        # Peripheral distribution - increases toward edge
        x = np.linspace(0, 1, n_bins)
        profiles.append(0.3 + 0.7 * x + np.random.randn(n_bins) * 0.05)

        # Perinuclear distribution - peaks early, decreases
        profiles.append(np.exp(-x * 2) + np.random.randn(n_bins) * 0.05)

        # Complex profile with multiple peaks
        profiles.append(0.5 + 0.3 * np.sin(x * 4 * np.pi) + 0.2 * np.cos(x * 2 * np.pi) + np.random.randn(n_bins) * 0.05)

        data_matrix = np.array(profiles)

        # Compare original vs vectorized
        results_orig = np.apply_along_axis(find_end_slope2, 1, data_matrix)
        results_vec = find_end_slope2_vectorized(data_matrix)

        np.testing.assert_array_equal(results_orig[:, 0], results_vec[:, 0])
        np.testing.assert_allclose(results_orig[:, 1], results_vec[:, 1], rtol=1e-10)

    def test_output_shape(self):
        """Test that output shape is correct."""
        n_rows = 50
        n_cols = 16
        data = np.random.randn(n_rows, n_cols)

        results = find_end_slope2_vectorized(data)

        assert results.shape == (n_rows, 2), f"Expected shape ({n_rows}, 2), got {results.shape}"

    def test_large_batch_performance(self):
        """Test performance with a large batch of profiles."""
        # This test mainly ensures the vectorized version can handle large batches
        # without errors (actual performance comparison would be done in a benchmark)
        n_rows = 10000
        n_cols = 16

        np.random.seed(42)
        data = np.random.randn(n_rows, n_cols)

        # Should complete without error
        results = find_end_slope2_vectorized(data)

        assert results.shape == (n_rows, 2)
        assert not np.any(np.isnan(results)), "Results contain NaN values"
        assert not np.any(np.isinf(results)), "Results contain infinite values"


@pytest.mark.benchmark
def test_performance_comparison():
    """Benchmark comparison between original and vectorized implementations.

    Run with: pytest -v -m benchmark test_vectorized_slope.py
    """
    import time

    n_rows = 1000
    n_cols = 16

    np.random.seed(42)
    data = np.random.randn(n_rows, n_cols)

    # Time original implementation
    start = time.time()
    results_orig = np.apply_along_axis(find_end_slope2, 1, data)
    time_orig = time.time() - start

    # Time vectorized implementation
    start = time.time()
    results_vec = find_end_slope2_vectorized(data)
    time_vec = time.time() - start

    # Verify results match
    np.testing.assert_array_equal(results_orig[:, 0], results_vec[:, 0])
    np.testing.assert_allclose(results_orig[:, 1], results_vec[:, 1], rtol=1e-10)

    speedup = time_orig / time_vec
    print(f"\nPerformance comparison ({n_rows} rows):")
    print(f"  Original (apply_along_axis): {time_orig:.4f}s")
    print(f"  Vectorized: {time_vec:.4f}s")
    print(f"  Speedup: {speedup:.1f}x")

    # Expect at least 5x speedup
    assert speedup > 5, f"Expected >5x speedup, got {speedup:.1f}x"
