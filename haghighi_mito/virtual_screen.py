"""Core virtual screen pipeline for mitochondrial morphology analysis.

This module contains the production pipeline for virtual screening across multiple
perturbation datasets (compounds, ORFs, CRISPR). For diagnostic and visualization
functions, see haghighi_mito.diagnostics.

Pipeline Overview
-----------------
1. Load per-site profiles and metadata
2. Calculate slopes from radial distribution patterns (bins 5-16)
3. Z-score normalize slopes and peak indices per plate
4. Aggregate via median across plates
5. Calculate statistical tests (Hotelling's T², Welch's t-test, Cohen's d)

Output Metrics
--------------
The pipeline generates these metrics per perturbation:

Statistical Tests:
- t_target_pattern: Hotelling's T² on full radial distribution (bins 5-16)
- t_orth: Hotelling's T² on orthogonal features (non-radial measurements)
- t_slope: Welch's t-test on z-scored slope values
- d_slope: Cohen's d effect size for slope

Direct Metrics:
- Count_Cells_avg: Mean cell count per perturbation
- last_peak_ind: Z-scored peak index (median across plates)
- slope: Z-scored slope from peak to end (median across plates)

Note: last_peak_ind and slope are z-score normalized per plate before aggregation,
matching the baseline methodology (notebook 2.0, lines 1251-1254).

Usage:
    pixi run haghighi-mito virtual-screen --dataset taorf
    pixi run haghighi-mito virtual-screen --dataset jump_orf --compare-baseline

    # For diagnostics/visualization:
    pixi run haghighi-mito compare-baseline-metrics --dataset taorf
"""

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger
from scipy.signal import savgol_filter

from haghighi_mito.config import (
    DATASET_INFO,
    EXTERNAL_DATA_DIR,
    PROCESSED_DATA_DIR,
)
from haghighi_mito.vectorized_slope import find_end_slope2_vectorized
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
    4. Z-score normalize slope and last_peak_ind per plate
    5. Aggregate per perturbation
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
    logger.info("Subtracting controls and calculating slopes (vectorized)...")

    # Build matrix of radial patterns (all rows at once)
    radial_patterns = per_site_df[target_columns].values  # shape: (n_rows, 12)

    # Build control matrix - map batch_plate to control values
    control_matrix = np.zeros_like(radial_patterns)
    batch_plates = per_site_df["batch_plate"].values

    for i, batch_plate in enumerate(batch_plates):
        if batch_plate in control_df_perplate.index:
            control_matrix[i] = control_df_perplate.loc[batch_plate].values

    # Subtract controls (vectorized)
    radial_patterns_corrected = radial_patterns - control_matrix

    # VECTORIZED SLOPE CALCULATION - processes all rows at once (~200x faster)
    slope_results = find_end_slope2_vectorized(radial_patterns_corrected)

    per_site_df["last_peak_ind"] = slope_results[:, 0]
    per_site_df["last_peak_loc"] = slope_results[:, 0]  # Alias for batch_plate_statistics compatibility
    per_site_df["slope"] = slope_results[:, 1]

    # Reset index to ensure clean indexing
    per_site_df = per_site_df.reset_index(drop=True)

    logger.info(f"Calculated slopes for {len(slope_results)} observations")

    # Step 3: Z-score normalize slope and last_peak_ind per plate
    # This matches notebook 2.0 line 1251-1254: standardize_per_catX normalizes per batch_plate
    logger.info("Z-score normalizing slope and last_peak_ind per plate...")

    def z_score_normalize(group):
        """Z-score normalize columns within a group."""
        result = group.copy()
        for col in ["slope", "last_peak_ind"]:
            mean_val = group[col].mean()
            std_val = group[col].std()
            if std_val > 0:  # Avoid division by zero
                result[col] = (group[col] - mean_val) / std_val
            else:
                result[col] = 0.0
        return result

    per_site_df = per_site_df.groupby("batch_plate", group_keys=False).apply(z_score_normalize)
    # Also update last_peak_loc alias to match
    per_site_df["last_peak_loc"] = per_site_df["last_peak_ind"]

    logger.info("Z-score normalization complete")

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


def calculate_statistical_tests(per_site_df, dataset: str):
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

    Returns
    -------
    pd.DataFrame
        Results with t-values for each perturbation
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

    logger.info(f"Calculated statistical tests for {len(results_list)} perturbations")

    return pd.DataFrame(results_list)


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
    output_dir = PROCESSED_DATA_DIR / "virtual_screen_module"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{dataset}_results_pattern_aug_070624.csv"
    results.to_csv(output_path, index=False)
    logger.info(f"\nSaved results to {output_path}")

    # Compare with baseline if requested
    if compare_baseline:
        from haghighi_mito.diagnostics import compare_with_baseline
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
