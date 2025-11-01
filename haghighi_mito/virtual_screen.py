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
    # Run virtual screen analysis
    pixi run haghighi-mito virtual-screen --dataset taorf

    # For diagnostic commands (compare with baseline, generate plots), see:
    pixi run haghighi-mito --help
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


def preprocess_metadata(dataset: str, mito_project_root: Path) -> pd.DataFrame:
    """Preprocess raw metadata for a dataset.

    This replicates the preprocessing logic from notebook 2.0-mh-virtual-screen.py
    (lines 132-262) to create standardized metadata with Batch, batch_plate, and
    ctrl_well columns.

    Args:
        dataset: Dataset name
        mito_project_root: Path to mito_project root directory

    Returns:
        Preprocessed metadata DataFrame
    """
    logger.info(f"Preprocessing metadata for {dataset}...")

    if dataset == "lincs":
        # Load raw metadata
        annot0 = pd.read_csv(mito_project_root / "workspace/metadata/lincs/DrugRepurposing_Metadata.csv")
        annot = pd.read_csv(mito_project_root / "workspace/metadata/LINCS_meta.csv")

        # Merge perturbation names
        annot = annot.merge(
            annot0[["Metadata_Plate", "Metadata_Well", "Metadata_pert_name"]],
            how="left",
            on=["Metadata_Plate", "Metadata_Well"],
        )

        # Add standard columns
        annot["Batch"] = "2016_04_01_a549_48hr_batch1_Mito_Project"
        annot["batch_plate"] = annot["Batch"] + "-" + annot["Metadata_Plate"]
        annot["ctrl_well"] = annot["Metadata_pert_type"].isin(["control"])

    elif dataset == "taorf":
        # Load raw metadata
        annot = pd.read_csv(mito_project_root / "workspace/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz")

        # Add standard columns
        annot["Batch"] = "2013_10_11_SIGMA2_Pilot"
        annot["batch_plate"] = annot["Batch"] + "-" + annot["Metadata_Plate"].astype(str)
        annot["ctrl_well"] = annot["Metadata_gene_name"].isin(["LacZ", "Luciferase"])
        annot["Metadata_pert_type"] = annot["Metadata_ASSAY_WELL_ROLE"]

        # Keep only necessary columns (as in notebook line 244-256)
        annot = annot[
            [
                "Metadata_Plate",
                "Metadata_Well",
                "Metadata_gene_name",
                "Metadata_pert_name",
                "Metadata_pert_type",
                "Metadata_broad_sample",
                "Metadata_moa",
                "batch_plate",
                "Batch",
                "ctrl_well",
            ]
        ]

    elif dataset == "CDRP":
        annot = pd.read_csv(mito_project_root / "workspace/metadata/CDRP_meta.csv")
        annot["Batch"] = "CDRP"
        annot["batch_plate"] = annot["Batch"] + "-" + annot["Metadata_Plate"].astype(str)
        annot["ctrl_well"] = annot["Metadata_pert_type"].isin(["control"])

    elif dataset == "jump_orf":
        plates = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/plate.csv.gz")
        wells = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/well.csv.gz")
        orf = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/orf.csv.gz")

        annot = wells.merge(plates, on="Metadata_PlateID", how="left")
        annot = annot.merge(orf, on="Metadata_JCP2022", how="left")
        annot["batch_plate"] = annot["Metadata_Batch"] + "-" + annot["Metadata_Plate"]
        annot["ctrl_well"] = annot["Metadata_control_type"].notna()
        annot["Metadata_pert_type"] = annot["Metadata_control_type"].fillna("trt")

    elif dataset == "jump_crispr":
        plates = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/plate.csv.gz")
        wells = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/well.csv.gz")
        crispr = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/crispr.csv.gz")

        annot = wells.merge(plates, on="Metadata_PlateID", how="left")
        annot = annot.merge(crispr, on="Metadata_JCP2022", how="left")
        annot["batch_plate"] = annot["Metadata_Batch"] + "-" + annot["Metadata_Plate"]
        annot["ctrl_well"] = annot["Metadata_control_type"].notna()
        annot["Metadata_pert_type"] = annot["Metadata_control_type"].fillna("trt")

    elif dataset == "jump_compound":
        plates = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/plate.csv.gz")
        wells = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/well.csv.gz")
        compound = pd.read_csv(mito_project_root / "workspace/metadata/JUMP/compound.csv.gz")

        annot = wells.merge(plates, on="Metadata_PlateID", how="left")
        annot = annot.merge(compound, on="Metadata_JCP2022", how="left")
        annot["batch_plate"] = annot["Metadata_Batch"] + "-" + annot["Metadata_Plate"]
        annot["ctrl_well"] = annot["Metadata_control_type"].notna()
        annot["Metadata_pert_type"] = annot["Metadata_control_type"].fillna("trt")

    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    # Ensure Metadata_Plate is string type (required for merging)
    if "Metadata_Plate" in annot.columns:
        annot["Metadata_Plate"] = annot["Metadata_Plate"].astype(str)

    logger.info(f"Preprocessed metadata: {len(annot)} rows")
    return annot


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

    # Preprocess metadata from raw files
    annot = preprocess_metadata(dataset, mito_project_root)

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


def calculate_metrics(per_site_df, annot, dataset: str):
    """Calculate core metrics from radial distribution data.

    Mimics notebook 2.0 logic:
    1. Calculate per-plate control means
    2. Subtract controls from target features
    3. Calculate slope per observation (vectorized)
    4. Z-score normalize slope and last_peak_ind per plate
    5. Aggregate per perturbation with all metadata columns

    Parameters
    ----------
    per_site_df : pd.DataFrame
        Per-site observations with all features
    annot : pd.DataFrame
        Annotation dataframe with all metadata columns
    dataset : str
        Dataset name

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (results_df with metadata + Count_Cells_avg + slope metrics, per_site_df with slopes added)
    """
    logger.info("Calculating core metrics...")

    pert_col = DATASET_INFO[dataset]["pert_col"]
    meta_cols = DATASET_INFO[dataset]["meta_cols"]

    # Define target columns (radial bins 5-16)
    target_columns = [f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16" for i in range(5, 17)]

    logger.info(f"Using {len(target_columns)} radial distribution bins (5-16)")

    # Step 1: Calculate per-plate control means
    logger.info("Calculating per-plate control means...")
    control_df_perplate = per_site_df.loc[per_site_df["ctrl_well"]].groupby("batch_plate")[target_columns].mean()

    logger.info(f"Found controls for {len(control_df_perplate)} plates")

    # Step 1b: Filter to only plates with controls (matches notebook lines 1230-1233)
    plates_with_controls = list(control_df_perplate.index)
    n_before = len(per_site_df)
    per_site_df = per_site_df[per_site_df["batch_plate"].isin(plates_with_controls)].reset_index(drop=True)
    logger.info(f"Filtered to plates with controls: {n_before} -> {len(per_site_df)} rows")

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

    # Use transform instead of apply to avoid FutureWarning about operating on grouping columns
    # transform is designed for per-group operations and automatically preserves all columns
    for col in ["slope", "last_peak_ind"]:
        per_site_df[col] = per_site_df.groupby("batch_plate")[col].transform(
            lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0.0
        )
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

    # START WITH ALL METADATA COLUMNS (matches notebook 2.0 line 1333)
    # This ensures we have Metadata_gene_name, Metadata_pert_name, Metadata_moa, etc.
    results = annot[meta_cols].drop_duplicates().reset_index(drop=True)
    logger.info(f"Starting with {len(results)} unique perturbations from metadata")

    # SINGLE-STAGE AGGREGATION
    # NOTE: This differs from the current notebook (which uses two-stage aggregation),
    # but empirically matches the July 2024 baseline better (r=0.90 vs 0.83 for slope).
    # The baseline appears to have been generated with single-stage aggregation.
    numeric_aggs = (
        pert_df.groupby(pert_col)
        .agg(
            {
                "Count_Cells": ["mean", "std", "count"],  # Average, variability, and count
                "last_peak_ind": "median",  # Median peak index
                "slope": ["median", "std"],  # Median slope + variability
                "batch_plate": "nunique",  # Number of unique plates
                "Metadata_Well": "nunique",  # Number of unique wells
            }
        )
        .reset_index()
    )

    # Flatten multi-level column names
    numeric_aggs.columns = [
        f"{col[0]}_{col[1]}" if col[1] else col[0] for col in numeric_aggs.columns
    ]

    # Rename to match expected output format
    numeric_aggs.rename(
        columns={
            "Count_Cells_mean": "Count_Cells_avg",
            "Count_Cells_std": "Count_Cells_std",
            "Count_Cells_count": "n_sites",  # Number of site observations
            "last_peak_ind_median": "last_peak_ind",  # Remove _median suffix
            "slope_median": "slope",
            "slope_std": "slope_std",
            "batch_plate_nunique": "n_plates",
            "Metadata_Well_nunique": "n_wells",
        },
        inplace=True,
    )

    # Merge numeric results into metadata results
    results = pd.merge(results, numeric_aggs, on=pert_col, how="left")

    logger.info(f"Aggregated results for {len(results)} perturbations")

    return results, per_site_df


def load_orthogonal_features(dataset: str):
    """Load list of orthogonal features for a dataset.

    Mirrors logic from notebook 2.0 lines 1005-1033:
    - LINCS uses hardcoded list
    - All other datasets use fibroblast_derived.csv
    """
    orth_features_dir = EXTERNAL_DATA_DIR / "mito_project/workspace/results/target_pattern_orth_features_lists"

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
    target_columns = [f"Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16" for i in range(5, 17)]

    # Load orthogonal features
    orth_features = load_orthogonal_features(dataset)

    # Prepare control data by plate
    logger.info("Preparing control data by plate...")
    control_df = per_site_df[per_site_df["ctrl_well"]].copy()
    # Dict comprehension is correct here (not dict(groupby)) - GroupBy objects aren't directly convertible to dict
    control_dfs_by_plate = {plate: group for plate, group in control_df.groupby("batch_plate")}
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
            logger.info(f"  Processed {i + 1}/{len(unique_perts)} perturbations")

        # Get data for this perturbation
        per_site_df_pert = pert_df[pert_df[pert_col] == pert].copy()

        # Calculate statistics using vectorized function
        batch_results = batch_plate_statistics(per_site_df_pert, control_dfs_by_plate, target_columns, orth_features)

        if batch_results is None:
            continue

        # Extract results
        tvals = batch_results["tvals"]  # shape (n_plates, 4)
        pvals = batch_results["pvals"]  # shape (n_plates, 6)
        plates = batch_results["plates"]

        # Find plate with median d_slope (index 3 of tvals) - using abs() sorting
        # NOTE: This differs from notebook (which uses percentile without abs())
        # but empirically matches baseline better (r=0.90 vs 0.85)
        median_plate_idx = np.argsort(np.abs(tvals[:, 3]))[len(tvals) // 2]
        median_plate_id = plates[median_plate_idx]

        # Get t-values from median plate
        t_target_pattern = tvals[median_plate_idx, 0]
        t_orth = tvals[median_plate_idx, 1]
        t_slope = tvals[median_plate_idx, 2]
        d_slope = tvals[median_plate_idx, 3]

        # Get p-values from median plate (matches notebook 2.0 lines 1401-1410)
        p_target_pattern = pvals[median_plate_idx, 0]
        p_orth = pvals[median_plate_idx, 1]
        p_slope = pvals[median_plate_idx, 2]
        p_slope_std = pvals[median_plate_idx, 3]
        p_pattern_std = pvals[median_plate_idx, 4]
        p_orth_std = pvals[median_plate_idx, 5]

        results_list.append(
            {
                pert_col: pert,
                "p_target_pattern": p_target_pattern,
                "p_orth": p_orth,
                "p_slope": p_slope,
                "p_slope_std": p_slope_std,
                "p_pattern_std": p_pattern_std,
                "p_orth_std": p_orth_std,
                "t_target_pattern": t_target_pattern,
                "t_orth": t_orth,
                "t_slope": t_slope,
                "d_slope": d_slope,
                "median_plate_id": median_plate_id,  # Track which plate was selected
            }
        )

    logger.info(f"Calculated statistical tests for {len(results_list)} perturbations")

    return pd.DataFrame(results_list)


def run_virtual_screen(dataset: str, calculate_stats: bool = True):
    """Run virtual screen analysis for specified dataset.

    Generates CSV with all metrics but does NOT compare with baseline.
    Use compare_with_baseline_csv() separately for comparison.

    Args:
        dataset: Dataset name (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)
        calculate_stats: Whether to calculate statistical tests (t-values, p-values)
    """
    if dataset not in DATASET_INFO:
        raise ValueError(f"Unknown dataset: {dataset}. Must be one of {list(DATASET_INFO.keys())}")

    logger.info(f"Starting virtual screen analysis for dataset: {dataset}")

    # Load data
    per_site_df, annot = load_dataset_data(dataset)

    # Calculate core metrics (metadata + Count_Cells_avg + slope metrics)
    results, per_site_df_with_slopes = calculate_metrics(per_site_df, annot, dataset)

    # Calculate statistical tests if requested
    if calculate_stats:
        stats_results = calculate_statistical_tests(per_site_df_with_slopes, dataset)

        # Merge with basic metrics
        pert_col = DATASET_INFO[dataset]["pert_col"]
        results = pd.merge(results, stats_results, on=pert_col, how="left")

    # Reorder columns to match baseline/notebook format:
    # Metadata cols, Count_Cells_avg, provenance metadata, p-values, t-values, last_peak_ind, slope
    meta_cols = DATASET_INFO[dataset]["meta_cols"]

    # Build column order
    column_order = meta_cols.copy()
    column_order.append("Count_Cells_avg")

    # Add provenance metadata (NEW - helps debug reproducibility issues)
    column_order.extend(
        [
            "n_sites",  # Number of site observations per perturbation
            "n_plates",  # Number of unique plates per perturbation
            "n_wells",  # Number of unique wells per perturbation
            "Count_Cells_std",  # Variability in cell counts across sites
            "slope_std",  # Variability in slope across sites
        ]
    )

    if calculate_stats:
        # Add p-value columns
        column_order.extend(
            [
                "p_target_pattern",
                "p_orth",
                "p_slope",
                "p_slope_std",
                "p_pattern_std",
                "p_orth_std",
            ]
        )
        # Add t-value columns
        column_order.extend(
            [
                "t_target_pattern",
                "t_orth",
                "t_slope",
                "d_slope",
                "median_plate_id",  # Which plate was selected for t-values
            ]
        )

    # Add slope columns
    column_order.extend(["last_peak_ind", "slope"])

    # Reorder (only include columns that exist in results)
    column_order = [col for col in column_order if col in results.columns]
    results = results[column_order]

    # Save results to symmetric location with baseline/regenerated structure
    output_dir = PROCESSED_DATA_DIR / "virtual_screen_module"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{dataset}_results_pattern_aug_070624.csv"
    results.to_csv(output_path, index=False)
    logger.info(f"\nSaved results to {output_path}")
