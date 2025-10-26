"""
Tests for vectorized statistical functions.

Verifies that vectorized implementations produce identical results to
the original loop-based implementations used in the virtual screen pipeline.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm, ttest_ind

from haghighi_mito.vectorized_stats import (
    TwoSampleT2Test_vectorized,
    batch_plate_statistics,
    cohens_d_vectorized,
    t_to_z_vectorized,
    ttest_ind_vectorized,
    z_to_p_vectorized,
)


# Original non-vectorized functions (from notebook) for comparison
def cohens_d_original(x, y):
    """Original Cohen's d implementation from notebook."""
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)
    return (np.mean(x) - np.mean(y)) / pooled_std


def t_to_z_original(t_stat, df):
    """Original t-to-z conversion from notebook."""
    return t_stat / np.sqrt(df / (df + t_stat**2))


def z_to_p_original(z):
    """Original z-to-p conversion from notebook."""
    return 2 * (1 - norm.cdf(abs(z)))


class TestCohensD:
    """Test vectorized Cohen's d calculation."""

    def test_single_comparison(self):
        """Test single comparison matches original."""
        np.random.seed(42)
        x = np.random.randn(20)
        y = np.random.randn(25)

        original = cohens_d_original(x, y)
        vectorized = cohens_d_vectorized([x], [y])[0]

        assert np.allclose(vectorized, original, rtol=1e-10)

    def test_multiple_comparisons(self):
        """Test multiple comparisons all match original."""
        np.random.seed(42)
        n_comparisons = 10

        x_vals = [np.random.randn(np.random.randint(10, 30)) for _ in range(n_comparisons)]
        y_vals = [np.random.randn(np.random.randint(10, 30)) for _ in range(n_comparisons)]

        original = [cohens_d_original(x, y) for x, y in zip(x_vals, y_vals, strict=False)]
        vectorized = cohens_d_vectorized(x_vals, y_vals)

        assert np.allclose(vectorized, original, rtol=1e-10)

    def test_equal_means(self):
        """Test Cohen's d is ~0 when means are equal."""
        np.random.seed(42)
        x = np.random.randn(50)
        y = np.random.randn(50)
        # Make means equal
        y = y - np.mean(y) + np.mean(x)

        result = cohens_d_vectorized([x], [y])[0]
        assert abs(result) < 0.1

    def test_large_effect(self):
        """Test Cohen's d detects large effects."""
        np.random.seed(42)
        x = np.random.randn(50)
        y = np.random.randn(50) + 2  # Large shift

        result = cohens_d_vectorized([x], [y])[0]
        assert abs(result) > 1.5  # Large effect size


class TestTTestVectorized:
    """Test vectorized t-test."""

    def test_single_ttest(self):
        """Test single t-test matches scipy."""
        np.random.seed(42)
        x = np.random.randn(20)
        y = np.random.randn(25) + 0.5

        original = ttest_ind(x, y, equal_var=False)
        t_stats, p_vals = ttest_ind_vectorized([x], [y], equal_var=False)

        assert np.allclose(t_stats[0], original.statistic, rtol=1e-10)
        assert np.allclose(p_vals[0], original.pvalue, rtol=1e-10)

    def test_multiple_ttests(self):
        """Test multiple t-tests match scipy."""
        np.random.seed(42)
        n_comparisons = 5

        x_vals = [np.random.randn(20) for _ in range(n_comparisons)]
        y_vals = [np.random.randn(25) + i * 0.3 for i in range(n_comparisons)]

        t_stats, p_vals = ttest_ind_vectorized(x_vals, y_vals, equal_var=False)

        for i in range(n_comparisons):
            original = ttest_ind(x_vals[i], y_vals[i], equal_var=False)
            assert np.allclose(t_stats[i], original.statistic, rtol=1e-10)
            assert np.allclose(p_vals[i], original.pvalue, rtol=1e-10)


class TestTToZConversion:
    """Test t-to-z conversion."""

    def test_vectorized_matches_original(self):
        """Test vectorized conversion matches original."""
        t_stats = np.array([1.5, 2.0, -1.2, 3.5, 0.0])
        dfs = np.array([10, 20, 30, 50, 100])

        original = np.array([t_to_z_original(t, df) for t, df in zip(t_stats, dfs, strict=False)])
        vectorized = t_to_z_vectorized(t_stats, dfs)

        assert np.allclose(vectorized, original, rtol=1e-10)

    def test_large_df_approaches_t(self):
        """Test that z approaches t for large degrees of freedom."""
        t_stat = 2.5
        large_df = 10000

        z = t_to_z_vectorized(np.array([t_stat]), np.array([large_df]))[0]

        # For very large df, z should be very close to t
        assert np.allclose(z, t_stat, rtol=0.01)


class TestZToP:
    """Test z-to-p conversion."""

    def test_vectorized_matches_original(self):
        """Test vectorized conversion matches original."""
        z_scores = np.array([0.0, 1.96, -1.96, 2.58, -2.58])

        original = np.array([z_to_p_original(z) for z in z_scores])
        vectorized = z_to_p_vectorized(z_scores)

        assert np.allclose(vectorized, original, rtol=1e-10)

    def test_standard_values(self):
        """Test known z-score to p-value conversions."""
        z_scores = np.array([1.96, 2.58])
        p_vals = z_to_p_vectorized(z_scores)

        # z=1.96 corresponds to p~0.05, z=2.58 to p~0.01
        assert np.allclose(p_vals[0], 0.05, rtol=0.01)
        assert np.allclose(p_vals[1], 0.01, rtol=0.02)  # Slightly relaxed tolerance


class TestTwoSampleT2Test:
    """Test vectorized Hotelling's T-squared test."""

    def test_identical_distributions(self):
        """Test that identical distributions give high p-value."""
        np.random.seed(42)
        X = np.random.randn(30, 5)
        Y = np.random.randn(30, 5)

        stats, p_vals, p_vals_std = TwoSampleT2Test_vectorized([X], [Y])

        # Should not be significant
        assert p_vals[0] > 0.05

    def test_different_distributions(self):
        """Test that different distributions give low p-value."""
        np.random.seed(42)
        X = np.random.randn(30, 5)
        Y = np.random.randn(30, 5) + 2  # Large shift

        stats, p_vals, p_vals_std = TwoSampleT2Test_vectorized([X], [Y])

        # Should be highly significant
        assert p_vals[0] < 0.001

    def test_multiple_tests(self):
        """Test multiple tests run correctly."""
        np.random.seed(42)
        n_tests = 3
        X_list = [np.random.randn(20, 4) for _ in range(n_tests)]
        Y_list = [np.random.randn(20, 4) + i for i in range(n_tests)]

        stats, p_vals, p_vals_std = TwoSampleT2Test_vectorized(X_list, Y_list)

        assert len(stats) == n_tests
        assert len(p_vals) == n_tests
        assert len(p_vals_std) == n_tests

        # Larger effect sizes should have smaller p-values
        assert p_vals[2] < p_vals[1] < p_vals[0]


class TestBatchPlateStatistics:
    """Test the main batch processing function."""

    def test_basic_functionality(self):
        """Test that batch_plate_statistics runs without errors."""
        np.random.seed(42)

        # Create synthetic data
        n_sites = 100
        plates = ["plate1", "plate2"]

        # Perturbation data
        per_site_df_pert = pd.DataFrame(
            {
                "batch_plate": np.repeat(plates, 50),
                "Count_Cells": np.random.randint(50, 200, n_sites),
                "last_peak_loc": np.random.randint(8, 16, n_sites),
                "slope": np.random.randn(n_sites) * 0.1 + 0.5,
            }
        )

        # Add target and orthogonal features
        for i in range(5, 17):
            per_site_df_pert[f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(n_sites) * 0.1

        for i in range(10):
            per_site_df_pert[f"Cells_feature_{i}"] = np.random.randn(n_sites) * 0.1

        # Control data
        control_dfs_by_plate = {}
        for plate in plates:
            control_dfs_by_plate[plate] = pd.DataFrame(
                {
                    "batch_plate": [plate] * 50,
                    "Count_Cells": np.random.randint(50, 200, 50),
                    "last_peak_loc": np.random.randint(8, 16, 50),
                    "slope": np.random.randn(50) * 0.1 + 0.3,
                }
            )

            for i in range(5, 17):
                control_dfs_by_plate[plate][f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(50) * 0.1

            for i in range(10):
                control_dfs_by_plate[plate][f"Cells_feature_{i}"] = np.random.randn(50) * 0.1

        target_columns = [f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16" for i in range(5, 17)]
        uncorr_feats = [f"Cells_feature_{i}" for i in range(10)]

        # Run batch statistics
        result = batch_plate_statistics(per_site_df_pert, control_dfs_by_plate, target_columns, uncorr_feats)

        assert result is not None
        assert "cell_counts" in result
        assert "pvals" in result
        assert "tvals" in result
        assert "peak_slope" in result
        assert "plates" in result

        assert len(result["cell_counts"]) == 2
        assert result["pvals"].shape == (2, 6)
        assert result["tvals"].shape == (2, 4)
        assert result["peak_slope"].shape == (2, 2)

    def test_single_plate(self):
        """Test with single plate."""
        np.random.seed(42)

        per_site_df_pert = pd.DataFrame(
            {
                "batch_plate": ["plate1"] * 30,
                "Count_Cells": np.random.randint(50, 200, 30),
                "last_peak_loc": np.random.randint(8, 16, 30),
                "slope": np.random.randn(30) * 0.1 + 0.5,
            }
        )

        for i in range(5, 17):
            per_site_df_pert[f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(30) * 0.1

        for i in range(5):
            per_site_df_pert[f"Cells_feature_{i}"] = np.random.randn(30) * 0.1

        control_dfs_by_plate = {
            "plate1": pd.DataFrame(
                {
                    "batch_plate": ["plate1"] * 30,
                    "Count_Cells": np.random.randint(50, 200, 30),
                    "last_peak_loc": np.random.randint(8, 16, 30),
                    "slope": np.random.randn(30) * 0.1 + 0.3,
                }
            )
        }

        for i in range(5, 17):
            control_dfs_by_plate["plate1"][f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(30) * 0.1

        for i in range(5):
            control_dfs_by_plate["plate1"][f"Cells_feature_{i}"] = np.random.randn(30) * 0.1

        target_columns = [f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16" for i in range(5, 17)]
        uncorr_feats = [f"Cells_feature_{i}" for i in range(5)]

        result = batch_plate_statistics(per_site_df_pert, control_dfs_by_plate, target_columns, uncorr_feats)

        assert result is not None
        assert len(result["plates"]) == 1
        assert result["pvals"].shape == (1, 6)


class TestEndToEnd:
    """End-to-end integration tests."""

    def test_full_pipeline_consistency(self):
        """Test that full pipeline produces consistent results."""
        np.random.seed(42)

        # Run twice with same seed
        results = []
        for _ in range(2):
            np.random.seed(42)

            per_site_df_pert = pd.DataFrame(
                {
                    "batch_plate": np.repeat(["plate1", "plate2"], 25),
                    "Count_Cells": np.random.randint(50, 200, 50),
                    "last_peak_loc": np.random.randint(8, 16, 50),
                    "slope": np.random.randn(50) * 0.1 + 0.5,
                }
            )

            for i in range(5, 17):
                per_site_df_pert[f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(50) * 0.1

            for i in range(5):
                per_site_df_pert[f"Cells_feature_{i}"] = np.random.randn(50) * 0.1

            control_dfs_by_plate = {}
            for plate in ["plate1", "plate2"]:
                control_dfs_by_plate[plate] = pd.DataFrame(
                    {
                        "batch_plate": [plate] * 25,
                        "Count_Cells": np.random.randint(50, 200, 25),
                        "last_peak_loc": np.random.randint(8, 16, 25),
                        "slope": np.random.randn(25) * 0.1 + 0.3,
                    }
                )

                for i in range(5, 17):
                    control_dfs_by_plate[plate][f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"] = np.random.randn(25) * 0.1

                for i in range(5):
                    control_dfs_by_plate[plate][f"Cells_feature_{i}"] = np.random.randn(25) * 0.1

            target_columns = [f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16" for i in range(5, 17)]
            uncorr_feats = [f"Cells_feature_{i}" for i in range(5)]

            result = batch_plate_statistics(per_site_df_pert, control_dfs_by_plate, target_columns, uncorr_feats)

            results.append(result)

        # Results should be identical
        assert np.allclose(results[0]["pvals"], results[1]["pvals"])
        assert np.allclose(results[0]["tvals"], results[1]["tvals"])
        assert np.allclose(results[0]["peak_slope"], results[1]["peak_slope"])
