"""Virtual screen analysis from scratch.

Starting point for rewriting the virtual screen pipeline with simple, verifiable logic.
Eventually will include orthogonal features, statistical testing, and all metrics.

Baseline Output Metrics
-----------------------
The July 2024 baseline CSVs contain these metrics per perturbation:

Statistical Tests (6 p-values, 4 t-values):
- p_target_pattern / t_target_pattern: Hotelling's T² test on entire radial distribution
  pattern (bins 5-16). Tests whether control-subtracted radial curves differ between
  treatment and control. Does NOT require peak detection or slope calculation.

- p_orth / t_orth: Hotelling's T² test on orthogonal features (non-radial measurements).
  Ensures perturbation doesn't affect unrelated cellular features.

- p_slope / t_slope: Welch's t-test on slope values (calculated from peak detection).
  Tests whether radial distribution slopes differ between treatment and control.

- d_slope: Cohen's d effect size for slope. Simple calculation:
  (mean_treatment - mean_control) / pooled_std

- p_slope_std, p_pattern_std, p_orth_std: Standardized p-values from z-score conversions.

Direct Metrics:
- Count_Cells_avg: Mean cell count per perturbation
- last_peak_ind: Index of last peak/valley in smoothed radial pattern (0-11 for bins 5-16)
- slope: Slope from last peak to end of radial distribution

Most Direct Baseline Comparison
--------------------------------
p_target_pattern bypasses peak detection entirely - it compares the full radial
distribution curves using multivariate statistics. This is the most fundamental test
and should match baseline if control subtraction is performed identically.

Usage:
    pixi run haghighi-mito virtual-screen --dataset taorf
    pixi run haghighi-mito virtual-screen --dataset jump_orf --compare-baseline
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from scipy.signal import savgol_filter

from haghighi_mito.config import (
    DATASET_INFO,
    EXTERNAL_DATA_DIR,
    PROCESSED_DATA_DIR,
)
from haghighi_mito.vectorized_stats import batch_plate_statistics


def find_end_slope2_simple(data):
    """Simplified version of find_end_slope2 from notebook 2.0.

    Finds the last peak/valley in the data and calculates slope from there to end.
    Returns (last_peak_index, slope).
    """
    # Smooth the data
    smoothed = savgol_filter(data, window_length=5, polyorder=3)

    # Find min/max indices
    min_max_indc = [np.argmax(smoothed), np.argmin(smoothed)]

    # Get indices that are not at the end (at least 2 positions before end)
    last_peak_ind0 = [i for i in min_max_indc if i < len(smoothed) - 2]

    if not last_peak_ind0:
        return 0, 0

    last_peak_ind = np.max(last_peak_ind0)

    # Calculate slope from last peak to average of last two points
    last_two_points_amplitude = (smoothed[-1] + smoothed[-2]) / 2
    slope = (last_two_points_amplitude - smoothed[last_peak_ind]) / (len(smoothed) - last_peak_ind - 1)

    return last_peak_ind, slope


def load_dataset_data(dataset: str):
    """Load per-site profiles and metadata for a given dataset.

    Args:
        dataset: Dataset name (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)

    Returns:
        Tuple of (per_site_df, annot)
    """
    logger.info(f"Loading {dataset} data...")

    # Paths
    mito_project_root = EXTERNAL_DATA_DIR / "mito_project"
    per_site_dir = mito_project_root / "workspace/per_site_aggregated_profiles_newpattern_2" / dataset
    metadata_path = mito_project_root / f"workspace/metadata/preprocessed/annot_{dataset}.csv"

    # Load metadata
    logger.info(f"Loading metadata from {metadata_path}")
    annot = pd.read_csv(metadata_path, dtype={"Metadata_Plate": str})

    # Load all per-site profile files in the directory
    logger.info(f"Loading per-site profiles from {per_site_dir}")
    per_site_files = list(per_site_dir.glob("*_site_agg_profiles.csv.gz"))

    if not per_site_files:
        raise FileNotFoundError(f"No per-site profile files found in {per_site_dir}")

    logger.info(f"Found {len(per_site_files)} per-site profile file(s)")

    per_site_dfs = []
    for per_site_file in per_site_files:
        logger.info(f"  Loading {per_site_file.name}")
        df = pd.read_csv(per_site_file)
        per_site_dfs.append(df)

    per_site_df = pd.concat(per_site_dfs, axis=0, ignore_index=True)

    # Convert Metadata_Plate to string first
    per_site_df["Metadata_Plate"] = per_site_df["Metadata_Plate"].astype(str)

    # Add batch_plate column
    per_site_df["batch_plate"] = per_site_df["Metadata_Batch"] + "-" + per_site_df["Metadata_Plate"]

    # Merge with metadata
    common_cols = list(set(annot.columns) & set(per_site_df.columns))

    # Use inner join for jump_crispr and jump_compound (as in notebook 2.0)
    merge_how = "inner" if dataset in ["jump_crispr", "jump_compound"] else "left"

    logger.info(f"Merging per-site profiles with metadata ({merge_how} join) on {len(common_cols)} columns")
    per_site_df = pd.merge(per_site_df, annot, how=merge_how, on=common_cols)

    # Filter out rows with null ctrl_well
    if "ctrl_well" in per_site_df.columns:
        n_before = len(per_site_df)
        per_site_df = per_site_df[~per_site_df["ctrl_well"].isnull()].reset_index(drop=True)
        logger.info(f"Filtered null ctrl_well: {n_before} -> {len(per_site_df)} rows")

    logger.info(f"Loaded {len(per_site_df)} per-site observations")
    logger.info(f"Unique perturbations: {per_site_df[DATASET_INFO[dataset]['pert_col']].nunique()}")

    return per_site_df, annot


def calculate_simple_metrics(per_site_df, dataset: str):
    """Calculate simple metrics from radial distribution data.

    Mimics notebook 2.0 logic but in minimal form:
    1. Calculate per-plate control means
    2. Subtract controls from target features
    3. Calculate slope per observation
    4. Aggregate per perturbation
    """
    logger.info("Calculating simple metrics...")

    pert_col = DATASET_INFO[dataset]["pert_col"]

    # Define target columns (radial bins 5-16)
    target_columns = [
        f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"
        for i in range(5, 17)
    ]

    logger.info(f"Using {len(target_columns)} radial distribution bins (5-16)")

    # Step 1: Calculate per-plate control means
    logger.info("Calculating per-plate control means...")
    control_df_perplate = (
        per_site_df.loc[per_site_df["ctrl_well"]]
        .groupby("batch_plate")[target_columns]
        .mean()
    )

    logger.info(f"Found controls for {len(control_df_perplate)} plates")

    # Step 2: Subtract controls per plate and calculate slope
    logger.info("Subtracting controls and calculating slopes...")

    slopes = []
    last_peak_inds = []

    for idx, row in per_site_df.iterrows():
        batch_plate = row["batch_plate"]

        # Get radial pattern for this observation
        radial_pattern = row[target_columns].values

        # Subtract control if available for this plate
        if batch_plate in control_df_perplate.index:
            control_values = control_df_perplate.loc[batch_plate].values
            radial_pattern_corrected = radial_pattern - control_values
        else:
            radial_pattern_corrected = radial_pattern

        # Calculate slope
        last_peak_ind, slope = find_end_slope2_simple(radial_pattern_corrected)
        last_peak_inds.append(last_peak_ind)
        slopes.append(slope)

    per_site_df["last_peak_ind"] = last_peak_inds
    per_site_df["last_peak_loc"] = last_peak_inds  # Alias for batch_plate_statistics compatibility
    per_site_df["slope"] = slopes

    # Reset index to ensure clean indexing
    per_site_df = per_site_df.reset_index(drop=True)

    logger.info(f"Calculated slopes for {len(slopes)} observations")

    # Step 4: Aggregate per perturbation (non-controls only)
    logger.info("Aggregating per perturbation...")

    # Create boolean masks
    is_not_control = per_site_df["ctrl_well"] == False  # noqa: E712
    has_pert_id = per_site_df[pert_col].notna()

    # Combine masks
    mask = is_not_control & has_pert_id

    # Filter using .loc
    pert_df = per_site_df.loc[mask, :].copy().reset_index(drop=True)

    logger.info(f"Filtered to {len(pert_df)} perturbation observations")
    logger.info(f"Unique perturbations: {pert_df[pert_col].nunique()}")

    results = pert_df.groupby(pert_col).agg({
        "Count_Cells": "mean",  # Average cell count
        "last_peak_ind": "median",  # Median peak index
        "slope": "median"  # Median slope
    }).reset_index()

    results.rename(columns={"Count_Cells": "Count_Cells_avg"}, inplace=True)

    logger.info(f"Aggregated results for {len(results)} perturbations")

    return results, per_site_df


def load_orthogonal_features(dataset: str):
    """Load list of orthogonal features for a dataset.

    Mirrors logic from notebook 2.0 lines 1005-1033:
    - LINCS uses hardcoded list
    - All other datasets use fibroblast_derived.csv
    """
    orth_features_dir = (
        EXTERNAL_DATA_DIR / "mito_project/workspace/results/target_pattern_orth_features_lists"
    )

    if dataset == "lincs":
        # Hardcoded for LINCS (from notebook 2.0, lines 1123-1132)
        orth_features = [
            "Nuclei_AreaShape_FormFactor",
            "Nuclei_AreaShape_Eccentricity",
            "Cells_AreaShape_Solidity",
            "Cells_Intensity_MaxIntensity_Mito",
            "Cells_AreaShape_Eccentricity",
            "Cytoplasm_AreaShape_MaxFeretDiameter",
            "Nuclei_Texture_AngularSecondMoment_DNA_8_45",
        ]
        logger.info(f"Using hardcoded orthogonal features for LINCS: {len(orth_features)} features")
    else:
        # Use fibroblast_derived.csv for all other datasets
        orth_features_path = orth_features_dir / "fibroblast_derived.csv"
        logger.info(f"Loading orthogonal features from {orth_features_path}")
        orth_features_df = pd.read_csv(orth_features_path)
        orth_features = orth_features_df["orth_fs"].tolist()
        logger.info(f"Loaded {len(orth_features)} orthogonal features")

    return orth_features


def calculate_statistical_tests(per_site_df, dataset: str, return_per_plate: bool = False):
    """Calculate statistical tests (t-values) for each perturbation.

    Uses the vectorized batch_plate_statistics function to compute:
    - t_target_pattern: Hotelling's T² for radial distribution
    - t_orth: Hotelling's T² for orthogonal features
    - t_slope: Welch's t-test for slope
    - d_slope: Cohen's d effect size for slope

    Parameters
    ----------
    per_site_df : pd.DataFrame
        Per-site dataframe with all observations
    dataset : str
        Dataset name
    return_per_plate : bool
        If True, return per-plate results instead of median-aggregated

    Returns
    -------
    pd.DataFrame
        Results with t-values for each perturbation (or per plate if return_per_plate=True)
    """
    logger.info("Calculating statistical tests...")

    pert_col = DATASET_INFO[dataset]["pert_col"]

    # Define target columns (radial bins 5-16)
    target_columns = [
        f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"
        for i in range(5, 17)
    ]

    # Load orthogonal features
    orth_features = load_orthogonal_features(dataset)

    # Prepare control data by plate
    logger.info("Preparing control data by plate...")
    control_df = per_site_df[per_site_df["ctrl_well"]].copy()
    control_dfs_by_plate = {
        plate: group for plate, group in control_df.groupby("batch_plate")
    }
    logger.info(f"Prepared controls for {len(control_dfs_by_plate)} plates")

    # Get unique perturbations (non-controls only)
    # Note: ctrl_well is dtype object, so use == False instead of ~ to avoid -1 indexing issue
    pert_df = per_site_df[per_site_df["ctrl_well"] == False].copy()  # noqa: E712
    unique_perts = pert_df[pert_col].dropna().unique()
    logger.info(f"Processing {len(unique_perts)} unique perturbations")

    # Process each perturbation
    results_list = []

    for i, pert in enumerate(unique_perts):
        if (i + 1) % 50 == 0:
            logger.info(f"  Processed {i+1}/{len(unique_perts)} perturbations")

        # Get data for this perturbation
        per_site_df_pert = pert_df[pert_df[pert_col] == pert].copy()

        # Calculate statistics using vectorized function
        batch_results = batch_plate_statistics(
            per_site_df_pert,
            control_dfs_by_plate,
            target_columns,
            orth_features
        )

        if batch_results is None:
            continue

        # Extract results
        tvals = batch_results['tvals']  # shape (n_plates, 4)
        plates = batch_results['plates']

        if return_per_plate:
            # Return one row per plate
            for pi in range(len(plates)):
                results_list.append({
                    pert_col: pert,
                    "batch_plate": plates[pi],
                    "t_target_pattern": tvals[pi, 0],
                    "t_orth": tvals[pi, 1],
                    "t_slope": tvals[pi, 2],
                    "d_slope": tvals[pi, 3],
                })
        else:
            # Find plate with median t_target_pattern (index 0)
            median_plate_idx = np.argsort(np.abs(tvals[:, 0]))[len(tvals) // 2]

            # Get t-values from median plate
            t_target_pattern = tvals[median_plate_idx, 0]
            t_orth = tvals[median_plate_idx, 1]
            t_slope = tvals[median_plate_idx, 2]
            d_slope = tvals[median_plate_idx, 3]

            results_list.append({
                pert_col: pert,
                "t_target_pattern": t_target_pattern,
                "t_orth": t_orth,
                "t_slope": t_slope,
                "d_slope": d_slope,
            })

    logger.info(f"Calculated statistical tests for {len(results_list)} {'plate-level' if return_per_plate else 'perturbation-level'} results")

    return pd.DataFrame(results_list)


def compare_with_baseline(results, dataset: str):
    """Compare calculated metrics with baseline CSV."""
    baseline_path = (
        EXTERNAL_DATA_DIR / "mito_project/workspace/results/virtual_screen_baseline"
        / f"{dataset}_results_pattern_aug_070624.csv"
    )

    logger.info(f"Loading baseline from {baseline_path}")
    baseline = pd.read_csv(baseline_path)

    # Merge on perturbation ID
    pert_col = DATASET_INFO[dataset]["pert_col"]

    # Determine which columns to compare based on what's in results
    baseline_cols = [pert_col, "Count_Cells_avg", "last_peak_ind", "slope"]

    # Add t-values if they're in results
    if "t_target_pattern" in results.columns:
        baseline_cols.extend(["t_target_pattern", "t_orth", "t_slope", "d_slope"])

    comparison = pd.merge(
        results,
        baseline[baseline_cols],
        on=pert_col,
        suffixes=("_new", "_baseline")
    )

    logger.info(f"Matched {len(comparison)} perturbations between new and baseline")

    # Calculate differences
    comparison["Count_Cells_diff"] = (
        comparison["Count_Cells_avg_new"] - comparison["Count_Cells_avg_baseline"]
    )
    comparison["Count_Cells_pct_diff"] = (
        100 * comparison["Count_Cells_diff"] / comparison["Count_Cells_avg_baseline"]
    )

    comparison["last_peak_ind_diff"] = (
        comparison["last_peak_ind_new"] - comparison["last_peak_ind_baseline"]
    )

    comparison["slope_diff"] = comparison["slope_new"] - comparison["slope_baseline"]
    comparison["slope_pct_diff"] = (
        100 * comparison["slope_diff"] / comparison["slope_baseline"].abs()
    )

    # Calculate t-value differences if present
    if "t_target_pattern_new" in comparison.columns:
        for col in ["t_target_pattern", "t_orth", "t_slope", "d_slope"]:
            comparison[f"{col}_diff"] = comparison[f"{col}_new"] - comparison[f"{col}_baseline"]
            comparison[f"{col}_pct_diff"] = (
                100 * comparison[f"{col}_diff"] / comparison[f"{col}_baseline"].abs()
            )

    # Summary statistics
    logger.info("\n" + "="*70)
    logger.info("COMPARISON WITH BASELINE")
    logger.info("="*70)

    logger.info(f"\nCount_Cells_avg:")
    logger.info(f"  Mean absolute diff: {comparison['Count_Cells_diff'].abs().mean():.4f}")
    logger.info(f"  Mean % diff: {comparison['Count_Cells_pct_diff'].abs().mean():.2f}%")
    logger.info(f"  Max % diff: {comparison['Count_Cells_pct_diff'].abs().max():.2f}%")
    logger.info(f"  Within 1%: {(comparison['Count_Cells_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    logger.info(f"\nlast_peak_ind:")
    logger.info(f"  Exact matches: {(comparison['last_peak_ind_diff'] == 0).sum()}/{len(comparison)}")
    logger.info(f"  Mean absolute diff: {comparison['last_peak_ind_diff'].abs().mean():.4f}")
    logger.info(f"  Max absolute diff: {comparison['last_peak_ind_diff'].abs().max():.0f}")

    logger.info(f"\nslope:")
    logger.info(f"  Mean absolute diff: {comparison['slope_diff'].abs().mean():.6f}")
    logger.info(f"  Mean % diff: {comparison['slope_pct_diff'].abs().mean():.2f}%")
    logger.info(f"  Max % diff: {comparison['slope_pct_diff'].abs().max():.2f}%")
    logger.info(f"  Within 10%: {(comparison['slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
    logger.info(f"  Within 1%: {(comparison['slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    # T-value statistics if present
    if "t_target_pattern_new" in comparison.columns:
        logger.info(f"\nt_target_pattern:")
        logger.info(f"  Mean absolute diff: {comparison['t_target_pattern_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_target_pattern_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_target_pattern_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_target_pattern_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info(f"\nt_orth:")
        logger.info(f"  Mean absolute diff: {comparison['t_orth_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_orth_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_orth_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_orth_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info(f"\nt_slope:")
        logger.info(f"  Mean absolute diff: {comparison['t_slope_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_slope_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info(f"\nd_slope:")
        logger.info(f"  Mean absolute diff: {comparison['d_slope_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['d_slope_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['d_slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['d_slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    logger.info("="*70)

    return comparison


def analyze_t_target_pattern_distribution(dataset: str):
    """Analyze the distribution of t_target_pattern values and their relationship to baseline.

    This function explores:
    1. Distribution statistics (mean, median, percentiles)
    2. Correlation between baseline and regenerated t_target_pattern
    3. Systematic transformations (linear, multiplicative)
    4. Relationship between t_target_pattern match quality and slope match quality

    Args:
        dataset: Dataset name
    """
    logger.info(f"Analyzing t_target_pattern distribution for {dataset}")

    # Load comparison CSV
    comparison_path = PROCESSED_DATA_DIR / "virtual_screen_simple" / f"{dataset}_baseline_comparison.csv"

    if not comparison_path.exists():
        logger.error(f"Comparison file not found: {comparison_path}")
        logger.error("Run virtual-screen with --compare-baseline first")
        return None

    comparison = pd.read_csv(comparison_path)
    logger.info(f"Loaded {len(comparison)} perturbations")

    # Check if t_target_pattern exists
    if "t_target_pattern_new" not in comparison.columns:
        logger.error("t_target_pattern not found in comparison data. Run with calculate_stats=True")
        return None

    logger.info(f"\n{'='*70}")
    logger.info("T_TARGET_PATTERN DISTRIBUTION ANALYSIS")
    logger.info(f"{'='*70}")

    # 1. Basic distribution statistics
    logger.info("\n1. DISTRIBUTION STATISTICS")
    logger.info("-" * 70)

    for suffix in ["new", "baseline"]:
        col = f"t_target_pattern_{suffix}"
        logger.info(f"\n{suffix.upper()} t_target_pattern:")
        logger.info(f"  Mean: {comparison[col].mean():.4f}")
        logger.info(f"  Median: {comparison[col].median():.4f}")
        logger.info(f"  Std: {comparison[col].std():.4f}")
        logger.info(f"  Min: {comparison[col].min():.4f}")
        logger.info(f"  Max: {comparison[col].max():.4f}")
        logger.info(f"  25th percentile: {comparison[col].quantile(0.25):.4f}")
        logger.info(f"  75th percentile: {comparison[col].quantile(0.75):.4f}")

    # 2. Correlation analysis
    logger.info("\n2. CORRELATION ANALYSIS")
    logger.info("-" * 70)

    from scipy.stats import pearsonr, spearmanr

    # Filter out NaN/inf values for correlation
    valid_mask = (
        np.isfinite(comparison["t_target_pattern_new"]) &
        np.isfinite(comparison["t_target_pattern_baseline"])
    )
    n_valid = valid_mask.sum()
    n_invalid = len(comparison) - n_valid

    if n_invalid > 0:
        logger.info(f"  Warning: {n_invalid} perturbations have NaN/inf values, using {n_valid} for correlation")

    if n_valid > 2:
        pearson_corr, pearson_p = pearsonr(
            comparison.loc[valid_mask, "t_target_pattern_new"],
            comparison.loc[valid_mask, "t_target_pattern_baseline"]
        )
        spearman_corr, spearman_p = spearmanr(
            comparison.loc[valid_mask, "t_target_pattern_new"],
            comparison.loc[valid_mask, "t_target_pattern_baseline"]
        )

        logger.info(f"  Pearson correlation: {pearson_corr:.4f} (p={pearson_p:.2e})")
        logger.info(f"  Spearman correlation: {spearman_corr:.4f} (p={spearman_p:.2e})")
    else:
        logger.info(f"  Insufficient valid data for correlation (n={n_valid})")

    # 3. Linear transformation analysis
    logger.info("\n3. LINEAR TRANSFORMATION ANALYSIS")
    logger.info("-" * 70)

    from scipy.stats import linregress

    if n_valid > 2:
        slope, intercept, r_value, p_value, std_err = linregress(
            comparison.loc[valid_mask, "t_target_pattern_baseline"],
            comparison.loc[valid_mask, "t_target_pattern_new"]
        )

        logger.info(f"  Linear fit: new = {slope:.4f} * baseline + {intercept:.4f}")
        logger.info(f"  R-squared: {r_value**2:.4f}")
        logger.info(f"  P-value: {p_value:.2e}")

        # Calculate how well linear transformation would work
        predicted = slope * comparison.loc[valid_mask, "t_target_pattern_baseline"] + intercept
        residuals = comparison.loc[valid_mask, "t_target_pattern_new"] - predicted
        logger.info(f"  Residual mean: {residuals.mean():.6f}")
        logger.info(f"  Residual std: {residuals.std():.6f}")
    else:
        logger.info(f"  Insufficient valid data for linear regression (n={n_valid})")

    # 4. Match quality comparison
    logger.info("\n4. MATCH QUALITY: t_target_pattern vs slope")
    logger.info("-" * 70)

    # Categorize by t_target_pattern match quality
    t_match_bins = [0, 10, 25, 50, 100, float('inf')]
    t_match_labels = ['<10%', '10-25%', '25-50%', '50-100%', '>100%']

    comparison['t_match_category'] = pd.cut(
        comparison['t_target_pattern_pct_diff'].abs(),
        bins=t_match_bins,
        labels=t_match_labels
    )

    for category in t_match_labels:
        mask = comparison['t_match_category'] == category
        n = mask.sum()
        if n == 0:
            continue

        slope_pct_mean = comparison.loc[mask, 'slope_pct_diff'].abs().mean()
        slope_within_10 = (comparison.loc[mask, 'slope_pct_diff'].abs() < 10).sum()

        logger.info(f"\n  {category} t_target_pattern match ({n} perturbations):")
        logger.info(f"    Mean slope % diff: {slope_pct_mean:.2f}%")
        logger.info(f"    Slope within 10%: {slope_within_10}/{n} ({100*slope_within_10/n:.1f}%)")

    # 5. Divergence patterns
    logger.info("\n5. DIVERGENCE PATTERNS")
    logger.info("-" * 70)

    # Find perturbations where t_target_pattern matches but slope diverges
    good_t_bad_slope = (
        (comparison['t_target_pattern_pct_diff'].abs() < 25) &
        (comparison['slope_pct_diff'].abs() > 50)
    )

    # Find perturbations where both match well
    both_good = (
        (comparison['t_target_pattern_pct_diff'].abs() < 25) &
        (comparison['slope_pct_diff'].abs() < 25)
    )

    # Find perturbations where both diverge
    both_bad = (
        (comparison['t_target_pattern_pct_diff'].abs() > 50) &
        (comparison['slope_pct_diff'].abs() > 50)
    )

    logger.info(f"  Good t_target_pattern (<25%), bad slope (>50%): {good_t_bad_slope.sum()}")
    logger.info(f"  Both good (<25%): {both_good.sum()}")
    logger.info(f"  Both bad (>50%): {both_bad.sum()}")

    # 6. Sign agreement
    logger.info("\n6. SIGN AGREEMENT (positive/negative)")
    logger.info("-" * 70)

    same_sign = (
        (comparison['t_target_pattern_new'] * comparison['t_target_pattern_baseline']) > 0
    )
    logger.info(f"  Same sign: {same_sign.sum()}/{len(comparison)} ({100*same_sign.mean():.1f}%)")

    same_slope_sign = (
        (comparison['slope_new'] * comparison['slope_baseline']) > 0
    )
    logger.info(f"  Same slope sign: {same_slope_sign.sum()}/{len(comparison)} ({100*same_slope_sign.mean():.1f}%)")

    # 7. Examples of perturbations where t_target_pattern matches but slope diverges
    logger.info("\n7. EXAMPLES: Good t_target_pattern, bad slope")
    logger.info("-" * 70)

    pert_col = DATASET_INFO[dataset]["pert_col"]

    # Get perturbations with good t_target_pattern (<25%) but bad slope (>50%)
    good_t_examples = comparison[good_t_bad_slope].nsmallest(10, "t_target_pattern_pct_diff")

    if len(good_t_examples) > 0:
        logger.info(f"\nTop 10 perturbations with best t_target_pattern match but poor slope match:")
        for idx, row in good_t_examples.iterrows():
            logger.info(f"\n  {row[pert_col]}:")
            logger.info(f"    t_target_pattern: new={row['t_target_pattern_new']:.4f}, baseline={row['t_target_pattern_baseline']:.4f}, diff={row['t_target_pattern_pct_diff']:.1f}%")
            logger.info(f"    slope: new={row['slope_new']:.6f}, baseline={row['slope_baseline']:.6f}, diff={row['slope_pct_diff']:.1f}%")
            logger.info(f"    last_peak_ind: new={row['last_peak_ind_new']:.1f}, baseline={row['last_peak_ind_baseline']:.1f}")
    else:
        logger.info("  No perturbations found with good t_target_pattern but bad slope")

    # 8. Create scatter plots
    logger.info("\n8. CREATING SCATTER PLOTS")
    logger.info("-" * 70)

    output_dir = PROCESSED_DATA_DIR / "figures" / "t_target_pattern_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create 2x2 subplot figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: t_target_pattern
    ax = axes[0, 0]
    if n_valid > 0:
        ax.scatter(
            comparison.loc[valid_mask, "t_target_pattern_baseline"],
            comparison.loc[valid_mask, "t_target_pattern_new"],
            alpha=0.5,
            s=20
        )
        # Add identity line
        min_val = min(comparison.loc[valid_mask, "t_target_pattern_baseline"].min(),
                     comparison.loc[valid_mask, "t_target_pattern_new"].min())
        max_val = max(comparison.loc[valid_mask, "t_target_pattern_baseline"].max(),
                     comparison.loc[valid_mask, "t_target_pattern_new"].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

        # Add linear fit line
        x_fit = np.linspace(min_val, max_val, 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, 'g--', alpha=0.7, label=f'y={slope:.2f}x+{intercept:.2f}')

        ax.set_xlabel('Baseline t_target_pattern', fontsize=11)
        ax.set_ylabel('Regenerated t_target_pattern', fontsize=11)
        ax.set_title(f't_target_pattern (r={pearson_corr:.3f}, R²={r_value**2:.3f})', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 2: slope
    ax = axes[0, 1]
    valid_slope_mask = (
        np.isfinite(comparison["slope_new"]) &
        np.isfinite(comparison["slope_baseline"])
    )
    if valid_slope_mask.sum() > 0:
        from scipy.stats import pearsonr as pr
        slope_corr, _ = pr(
            comparison.loc[valid_slope_mask, "slope_baseline"],
            comparison.loc[valid_slope_mask, "slope_new"]
        )

        ax.scatter(
            comparison.loc[valid_slope_mask, "slope_baseline"],
            comparison.loc[valid_slope_mask, "slope_new"],
            alpha=0.5,
            s=20
        )
        # Add identity line
        min_val = min(comparison.loc[valid_slope_mask, "slope_baseline"].min(),
                     comparison.loc[valid_slope_mask, "slope_new"].min())
        max_val = max(comparison.loc[valid_slope_mask, "slope_baseline"].max(),
                     comparison.loc[valid_slope_mask, "slope_new"].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

        ax.set_xlabel('Baseline slope', fontsize=11)
        ax.set_ylabel('Regenerated slope', fontsize=11)
        ax.set_title(f'slope (r={slope_corr:.3f})', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 3: t_orth
    ax = axes[1, 0]
    if "t_orth_new" in comparison.columns:
        valid_torth_mask = (
            np.isfinite(comparison["t_orth_new"]) &
            np.isfinite(comparison["t_orth_baseline"])
        )
        if valid_torth_mask.sum() > 0:
            torth_corr, _ = pr(
                comparison.loc[valid_torth_mask, "t_orth_baseline"],
                comparison.loc[valid_torth_mask, "t_orth_new"]
            )

            ax.scatter(
                comparison.loc[valid_torth_mask, "t_orth_baseline"],
                comparison.loc[valid_torth_mask, "t_orth_new"],
                alpha=0.5,
                s=20
            )
            # Add identity line
            min_val = min(comparison.loc[valid_torth_mask, "t_orth_baseline"].min(),
                         comparison.loc[valid_torth_mask, "t_orth_new"].min())
            max_val = max(comparison.loc[valid_torth_mask, "t_orth_baseline"].max(),
                         comparison.loc[valid_torth_mask, "t_orth_new"].max())
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

            ax.set_xlabel('Baseline t_orth', fontsize=11)
            ax.set_ylabel('Regenerated t_orth', fontsize=11)
            ax.set_title(f't_orth (r={torth_corr:.3f})', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)

    # Plot 4: t_slope
    ax = axes[1, 1]
    if "t_slope_new" in comparison.columns:
        valid_tslope_mask = (
            np.isfinite(comparison["t_slope_new"]) &
            np.isfinite(comparison["t_slope_baseline"])
        )
        if valid_tslope_mask.sum() > 0:
            tslope_corr, _ = pr(
                comparison.loc[valid_tslope_mask, "t_slope_baseline"],
                comparison.loc[valid_tslope_mask, "t_slope_new"]
            )

            ax.scatter(
                comparison.loc[valid_tslope_mask, "t_slope_baseline"],
                comparison.loc[valid_tslope_mask, "t_slope_new"],
                alpha=0.5,
                s=20
            )
            # Add identity line
            min_val = min(comparison.loc[valid_tslope_mask, "t_slope_baseline"].min(),
                         comparison.loc[valid_tslope_mask, "t_slope_new"].min())
            max_val = max(comparison.loc[valid_tslope_mask, "t_slope_baseline"].max(),
                         comparison.loc[valid_tslope_mask, "t_slope_new"].max())
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')

            ax.set_xlabel('Baseline t_slope', fontsize=11)
            ax.set_ylabel('Regenerated t_slope', fontsize=11)
            ax.set_title(f't_slope (r={tslope_corr:.3f})', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)

    fig.suptitle(f'{dataset.upper()}: Baseline vs Regenerated Metrics\n(n={len(comparison)} perturbations)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    plot_path = output_dir / f"{dataset}_baseline_vs_regenerated.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    logger.info(f"  Saved scatter plots to {plot_path}")
    plt.close()

    logger.info(f"\n{'='*70}")

    return comparison


def compare_per_plate_results(dataset: str):
    """Compare per-plate statistical results with baseline.

    For perturbations appearing on only ONE plate, the per-plate result
    must equal the aggregated result. This tests if our core per-plate
    calculations match the baseline without aggregation confounds.

    Args:
        dataset: Dataset name
    """
    logger.info(f"Comparing per-plate results for {dataset}")

    # Load data
    per_site_df, _ = load_dataset_data(dataset)

    # Calculate slopes first (needed by statistical tests)
    _, per_site_df_with_slopes = calculate_simple_metrics(per_site_df, dataset)

    # Calculate per-plate statistics
    logger.info("Calculating per-plate statistics (this may take a few minutes)...")
    per_plate_results = calculate_statistical_tests(per_site_df_with_slopes, dataset, return_per_plate=True)

    # Count plates per perturbation
    pert_col = DATASET_INFO[dataset]["pert_col"]
    plates_per_pert = per_plate_results.groupby(pert_col).size()
    single_plate_perts = plates_per_pert[plates_per_pert == 1].index.tolist()

    logger.info(f"Found {len(single_plate_perts)} perturbations appearing on exactly 1 plate")
    logger.info(f"Total perturbations: {len(plates_per_pert)}")

    # Load baseline
    baseline_path = (
        EXTERNAL_DATA_DIR / "mito_project/workspace/results/virtual_screen_baseline"
        / f"{dataset}_results_pattern_aug_070624.csv"
    )
    baseline = pd.read_csv(baseline_path)

    logger.info(f"\n{'='*70}")
    logger.info("PER-PLATE vs BASELINE COMPARISON")
    logger.info(f"{'='*70}")
    logger.info(f"\nFor each baseline value, checking if it matches ANY per-plate value")
    logger.info(f"(If baseline used different plate selection, values might match exactly)\n")

    # For each perturbation, check if baseline value matches any per-plate value
    matches = {"t_target_pattern": [], "t_orth": [], "t_slope": [], "d_slope": []}
    all_perts = baseline[pert_col].unique()

    for pert in all_perts:
        baseline_row = baseline[baseline[pert_col] == pert].iloc[0]
        pert_plates = per_plate_results[per_plate_results[pert_col] == pert]

        for metric in matches.keys():
            baseline_val = baseline_row[metric]
            # Check if baseline value matches any per-plate value within 0.01%
            match_found = any(abs(plate_val - baseline_val) / abs(baseline_val) < 0.0001
                             for plate_val in pert_plates[metric].values
                             if not pd.isna(plate_val) and baseline_val != 0)
            matches[metric].append(match_found)

    logger.info("Results: How many perturbations have baseline value matching ANY per-plate value?\n")

    for metric in ["t_target_pattern", "t_orth", "t_slope", "d_slope"]:
        n_matches = sum(matches[metric])
        pct_matches = 100 * n_matches / len(all_perts)
        logger.info(f"{metric}:")
        logger.info(f"  Exact matches (within 0.01%): {n_matches}/{len(all_perts)} ({pct_matches:.1f}%)")

    # Save per-plate results for inspection
    output_dir = PROCESSED_DATA_DIR / "virtual_screen_simple"
    output_path = output_dir / f"{dataset}_per_plate_results.csv"
    per_plate_results.to_csv(output_path, index=False)
    logger.info(f"\nSaved per-plate results to {output_path}")

    logger.info(f"\n{'='*70}")

    return per_plate_results


def run_virtual_screen(dataset: str, compare_baseline: bool = True, calculate_stats: bool = True):
    """Run virtual screen analysis for specified dataset.

    Args:
        dataset: Dataset name (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)
        compare_baseline: Whether to compare with baseline CSV and save comparison
        calculate_stats: Whether to calculate statistical tests (t-values)
    """
    if dataset not in DATASET_INFO:
        raise ValueError(f"Unknown dataset: {dataset}. Must be one of {list(DATASET_INFO.keys())}")

    logger.info(f"Starting virtual screen analysis for dataset: {dataset}")

    # Load data
    per_site_df, annot = load_dataset_data(dataset)

    # Calculate basic metrics
    results, per_site_df_with_slopes = calculate_simple_metrics(per_site_df, dataset)

    # Calculate statistical tests if requested
    if calculate_stats:
        stats_results = calculate_statistical_tests(per_site_df_with_slopes, dataset)

        # Merge with basic metrics
        pert_col = DATASET_INFO[dataset]["pert_col"]
        results = pd.merge(results, stats_results, on=pert_col, how="left")

    # Save results to symmetric location with baseline/regenerated structure
    output_dir = PROCESSED_DATA_DIR / "virtual_screen_simple"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{dataset}_results_pattern_aug_070624.csv"
    results.to_csv(output_path, index=False)
    logger.info(f"\nSaved results to {output_path}")

    # Compare with baseline if requested
    if compare_baseline:
        comparison = compare_with_baseline(results, dataset)

        # Save comparison
        comparison_path = output_dir / f"{dataset}_baseline_comparison.csv"
        comparison.to_csv(comparison_path, index=False)
        logger.info(f"\nSaved comparison to {comparison_path}")

        # Show some examples of large differences
        if calculate_stats and "t_target_pattern_pct_diff" in comparison.columns:
            logger.info("\nTop 5 largest t_target_pattern % differences:")
            pert_col = DATASET_INFO[dataset]["pert_col"]
            top_diffs = comparison.nlargest(5, "t_target_pattern_pct_diff")[
                [pert_col, "t_target_pattern_new", "t_target_pattern_baseline", "t_target_pattern_pct_diff"]
            ]
            print(top_diffs.to_string(index=False))
        else:
            logger.info("\nTop 5 largest slope % differences:")
            pert_col = DATASET_INFO[dataset]["pert_col"]
            top_diffs = comparison.nlargest(5, "slope_pct_diff")[
                [pert_col, "slope_new", "slope_baseline", "slope_pct_diff"]
            ]
            print(top_diffs.to_string(index=False))
