"""Patient fibroblast data loading and preprocessing.

Shared data loading pipeline used by notebooks 1.0 and 3.0. Reads single-cell
SQLite databases and pickle files, performs QC filtering, and returns
standardized DataFrames ready for downstream analysis.

The heavy path (SQLite + pickle) produces the same data that is available
pre-aggregated as aggregated_profiles_fibroblast.csv on S3.
"""

import os

import numpy as np
import pandas as pd
import sklearn.preprocessing as sp
from loguru import logger
from singlecell.preprocess.extract_cpfeature_names import extract_cpfeature_names
from singlecell.preprocess.find_highly_correlated_features import find_correlation
from singlecell.preprocess.handle_nans import handle_nans
from singlecell.read.read_single_cell_sql import readSingleCellData_sqlalch

from haghighi_mito.config import (
    AGGREGATED_PROFILES_PATH,
    FIBROBLAST_DATA_DIR,
    PATIENT_LABELS_PATH,
    SQLITE_DATA_DIR,
)


def load_patient_labels() -> pd.DataFrame:
    """Load and standardize patient disease labels.

    Returns DataFrame with columns: subject (str), D, D1, etc.
    """
    disease_labels = pd.read_excel(PATIENT_LABELS_PATH)
    disease_labels = disease_labels.rename(columns={"ID": "subject"})
    disease_labels["subject"] = disease_labels["subject"].astype(str)
    return disease_labels


def load_single_cell_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load single-cell data from SQLite databases and pickle files.

    Returns:
        (df_raw, df_allfeatures): Raw SQLite data merged with labels,
        and the allFeatures pickle DataFrame.
    """
    logger.info("Loading single-cell data from SQLite databases...")

    # Read all SQLite databases
    subject_list = sorted(os.listdir(SQLITE_DATA_DIR))
    # Filter to actual directories (skip .DS_Store etc.)
    subject_list = [s for s in subject_list if (SQLITE_DATA_DIR / s).is_dir()]
    logger.info(f"Found {len(subject_list)} subject directories")

    sc_df_ls = []
    compartments = ["Cells", "Cytoplasm", "Nuclei"]
    for si in subject_list:
        file_name = str(SQLITE_DATA_DIR / si / f"{si}.sqlite")
        if os.path.exists(file_name):
            sc_df = readSingleCellData_sqlalch(file_name, compartments)
            sc_df["subject"] = si
            sc_df_ls.append(sc_df)

    df_new = pd.concat(sc_df_ls, ignore_index=True)
    logger.info(f"Loaded {len(df_new)} cells from {len(sc_df_ls)} subjects")

    # Load pickle with labels and all features
    logger.info("Loading pickle files...")
    pkl_path = FIBROBLAST_DATA_DIR / "single_cell_with_annot_allFeatures.pkl"
    df_allfeatures = pd.read_pickle(pkl_path, compression="infer")
    df_allfeatures = df_allfeatures.interpolate()
    logger.info(f"allFeatures pickle: {df_allfeatures.shape}")

    # Merge: keep only subjects present in both
    common_subjects = list(set(df_new["subject"].unique()) & set(df_allfeatures["subject"].unique()))
    df_raw = df_new[df_new["subject"].isin(common_subjects)].copy()

    # Add labels from pickle
    df_raw = pd.merge(df_raw, df_allfeatures[["subject", "label"]].drop_duplicates(), on="subject", how="left")

    # Merge with disease labels
    disease_labels = load_patient_labels()
    df_raw = pd.merge(df_raw, disease_labels, on=["subject"], how="left")

    # Subject corrections
    df_raw["subject"] = df_raw["subject"].replace(["370E", "370F", "370H"], "370")
    df_raw.loc[df_raw["subject"] == "272", "label"] = "Control"
    df_raw.loc[df_raw["subject"] == "MCL004", "label"] = "SZA"

    logger.info(f"Merged data: {len(df_raw)} cells, {df_raw['subject'].nunique()} subjects")
    return df_raw, df_allfeatures


def preprocess_single_cells(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Apply QC filtering and feature selection to single-cell data.

    Steps:
    1. Extract CellProfiler feature names
    2. Remove NaN-heavy and low-variance features (handle_nans)
    3. Remove border cells (borderLength=200)
    4. Remove cells with bad cell:nucleus ratios
    5. Remove intensity artifact cells
    6. Second handle_nans pass
    7. Standardize features

    Returns:
        (df_raw, df_scaled, cp_features_analysis): Raw filtered DataFrame, standardized DataFrame,
        and list of analysis features. Both DataFrames have identical rows after QC filtering;
        df_raw preserves original values while df_scaled has z-scored features.
    """
    logger.info("Preprocessing single-cell data...")

    # Extract feature names
    cp_features, cp_features_analysis_0 = extract_cpfeature_names(df)
    df, cp_features_analysis = handle_nans(df, cp_features_analysis_0, thrsh_null_ratio=0.05, thrsh_std=0.001, fill_na_method="drop-rows")
    logger.info(f"After handle_nans: {df.shape}")

    # Remove border cells
    border_length = 200
    im_width = df["Width_Mito"].values[0]
    im_height = df["Height_Mito"].values[0]
    df = df.loc[
        ~(
            (df["Nuclei_Location_Center_X"] > (im_width - border_length))
            | (df["Nuclei_Location_Center_X"] < border_length)
            | (df["Nuclei_Location_Center_Y"] > (im_height - border_length))
            | (df["Nuclei_Location_Center_Y"] < border_length)
        ),
        :,
    ].reset_index(drop=True)
    logger.info(f"After border cell removal: {df.shape}")

    # Remove cells with bad segmentation (cell:nucleus ratio filters)
    df["Cells2Nuclei_MajorAxisLengthRatio"] = df["Cells_AreaShape_MajorAxisLength"] / df["Nuclei_AreaShape_MajorAxisLength"]
    df["Cells2Nuclei_AreaShapeRatio"] = df["Cells_AreaShape_Area"] / df["Nuclei_AreaShape_Area"]
    df = df[df["Cells2Nuclei_MajorAxisLengthRatio"] > 2].reset_index(drop=True)
    df = df[df["Cells2Nuclei_AreaShapeRatio"] > 5].reset_index(drop=True)
    logger.info(f"After cell seg ratio filter: {df.shape}")

    # Remove intensity artifact cells
    df = df[(df["Cells_Intensity_MeanIntensity_Actin"] < 0.5) & (df["Nuclei_Intensity_MeanIntensity_Actin"] < 0.55)].reset_index(drop=True)
    logger.info(f"After intensity artifact removal: {df.shape}")

    # Second handle_nans pass
    df, cp_features_analysis = handle_nans(df, cp_features_analysis, thrsh_null_ratio=0.05, thrsh_std=0.001, fill_na_method="drop-rows")

    # Standardize
    scaler = sp.StandardScaler()
    scaler.fit(df[cp_features_analysis])
    df["label"] = df["label"].replace("MDD or Dep", "MDD")
    df_scaled = df.copy()
    df_scaled[cp_features_analysis] = scaler.transform(df[cp_features_analysis])

    logger.info(f"Final processed data: {df_scaled.shape}, {len(cp_features_analysis)} features")
    return df, df_scaled, cp_features_analysis


def prepare_psychosis_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Add composite 'psychosis' group (BP + SZ + SZA) to DataFrame.

    Args:
        df: DataFrame with 'label' column containing patient categories.

    Returns:
        DataFrame with psychosis rows appended.
    """
    data = df.groupby(["label", "subject"]).mean(numeric_only=True).reset_index()
    psych = data[data["label"].isin(["BP", "SZ", "SZA"])].copy()
    psych["label"] = "psychosis"
    return pd.concat([data, psych], ignore_index=True)


def load_aggregated_profiles() -> pd.DataFrame:
    """Load pre-aggregated patient profiles with label corrections.

    This is the lightweight alternative to load_single_cell_data() +
    preprocess_single_cells(). Uses the pre-computed CSV from S3.
    """
    df = pd.read_csv(AGGREGATED_PROFILES_PATH)
    df.loc[df["subject"].astype(str) == "272", "label"] = "Control"
    df.loc[df["subject"].astype(str) == "MCL004", "label"] = "SZA"
    df["label"] = df["label"].replace("MDD or Dep", "MDD")
    return df
