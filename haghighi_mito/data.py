"""Data loading and database utilities."""

from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from loguru import logger

from haghighi_mito.config import DATASET_INFO, PROCESSED_DATA_DIR, PROCESSED_TABLES_DIR


def create_screen_database(
    output_path: Path = PROCESSED_DATA_DIR / "screen_results.duckdb",
    tables_dir: Path | None = None,
    overwrite: bool = False,
    use_parquet: bool = False,
    parquet_dir: Path | None = None,
) -> Path:
    """Combine all screen result Excel/Parquet files into a single DuckDB database.

    Only relevant columns are kept for each dataset to reduce database size.

    Args:
        output_path: Path to output DuckDB database file
        tables_dir: Directory containing Excel files. If None, uses PROCESSED_TABLES_DIR from config
        overwrite: If True, recreate database even if it exists
        use_parquet: If True, read from Parquet files instead of Excel
        parquet_dir: Directory containing Parquet files (required if use_parquet=True)

    Returns:
        Path to created database file
    """
    if output_path.exists() and not overwrite:
        logger.info(f"Database already exists at {output_path}. Use overwrite=True to recreate.")
        return output_path

    # Delete existing database if overwrite=True
    if output_path.exists():
        output_path.unlink()
        logger.info(f"Deleted existing database at {output_path}")

    # Validate parameters
    if use_parquet and parquet_dir is None:
        raise ValueError("parquet_dir must be provided when use_parquet=True")

    # Use provided tables_dir or default from config
    if tables_dir is None:
        tables_dir = PROCESSED_TABLES_DIR
    else:
        tables_dir = Path(tables_dir)

    # Set parquet_dir if using parquet
    if use_parquet:
        parquet_dir = Path(parquet_dir)

    # Common columns present in all datasets
    common_cols = ["d_slope", "p_slope_std", "p_orth_std", "t_orth", "Count_Cells_avg"]

    # Dataset-specific columns to keep
    dataset_cols = {
        "CDRP": common_cols
        + [
            "Metadata_Sample_Dose",
            "Metadata_broad_sample",
            "Metadata_pert_id",
            "Metadata_mmoles_per_liter2",
            "Metadata_moa",
        ],
        "JUMP_Compound": common_cols
        + ["Metadata_JCP2022", "Metadata_InChIKey", "Metadata_InChI"],
        "JUMP_CRISPR": common_cols
        + ["Metadata_JCP2022", "Metadata_Symbol", "Metadata_NCBI_Gene_ID"],
        "JUMP_ORF": common_cols
        + ["Metadata_JCP2022", "Metadata_Symbol", "Metadata_broad_sample"],
        "LINCS": common_cols
        + [
            "Metadata_pert_id_dose",
            "Metadata_pert_name",
            "Metadata_broad_sample",
            "Metadata_moa",
            "Metadata_target",
            "Metadata_InChIKey14",
        ],
        "TA_ORF": common_cols
        + ["Metadata_broad_sample", "Metadata_gene_name", "Metadata_pert_name", "Metadata_moa"],
    }

    # Mapping from dataset names to Excel sheet names (unfiltered sheets)
    sheet_names = {
        "CDRP": "CDRP",
        "JUMP_Compound": "jump_compound",
        "JUMP_CRISPR": "jump_crispr",
        "JUMP_ORF": "jump_orf",
        "LINCS": "lincs",
        "TA_ORF": "taorf",
    }

    # Screen result files
    if use_parquet:
        files = {
            "CDRP": parquet_dir / "CDRP_unfiltered.parquet",
            "JUMP_Compound": parquet_dir / "jump_compound_unfiltered.parquet",
            "JUMP_CRISPR": parquet_dir / "jump_crispr_unfiltered.parquet",
            "JUMP_ORF": parquet_dir / "jump_orf_unfiltered.parquet",
            "LINCS": parquet_dir / "lincs_unfiltered.parquet",
            "TA_ORF": parquet_dir / "taorf_unfiltered.parquet",
        }
    else:
        files = {
            "CDRP": tables_dir / "CDRP_screen_results.xlsx",
            "JUMP_Compound": tables_dir / "jump_compound_screen_results.xlsx",
            "JUMP_CRISPR": tables_dir / "jump_crispr_screen_results.xlsx",
            "JUMP_ORF": tables_dir / "jump_orf_screen_results.xlsx",
            "LINCS": tables_dir / "lincs_screen_results.xlsx",
            "TA_ORF": tables_dir / "taorf_screen_results.xlsx",
        }

    # Mapping for unified perturbation columns
    def add_perturbation_columns(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
        """Add Metadata_pert_type and Metadata_pert_id columns."""
        if dataset_name == "CDRP":
            df["Metadata_pert_type"] = "compound"
            df["Metadata_pert_id"] = df["Metadata_pert_id"]
        elif dataset_name == "JUMP_Compound":
            df["Metadata_pert_type"] = "compound"
            df["Metadata_pert_id"] = df["Metadata_InChIKey"]
        elif dataset_name == "JUMP_CRISPR":
            df["Metadata_pert_type"] = "gene_crispr"
            df["Metadata_pert_id"] = df["Metadata_Symbol"]
        elif dataset_name == "JUMP_ORF":
            df["Metadata_pert_type"] = "gene_orf"
            df["Metadata_pert_id"] = df["Metadata_Symbol"]
        elif dataset_name == "LINCS":
            df["Metadata_pert_type"] = "compound"
            df["Metadata_pert_id"] = df["Metadata_pert_name"]
        elif dataset_name == "TA_ORF":
            df["Metadata_pert_type"] = "gene_orf"
            df["Metadata_pert_id"] = df["Metadata_gene_name"]
        return df

    # Load and combine all datasets
    dfs = []
    for dataset_name, file_path in files.items():
        if use_parquet:
            logger.info(f"Loading {dataset_name} from {file_path.name}")
            df = pd.read_parquet(file_path)
        else:
            sheet_name = sheet_names[dataset_name]
            logger.info(f"Loading {dataset_name} from {file_path.name}, sheet '{sheet_name}'")
            df = pd.read_excel(file_path, sheet_name=sheet_name)

        # Select only relevant columns that exist in the dataframe
        cols_to_keep = [col for col in dataset_cols[dataset_name] if col in df.columns]
        df = df[cols_to_keep]

        # Add dataset source identifier
        df.insert(0, "Metadata_dataset", dataset_name)

        # Add unified perturbation columns
        df = add_perturbation_columns(df, dataset_name)

        logger.info(f"  Kept {len(cols_to_keep)} columns, {len(df):,} rows")
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined {len(combined):,} rows from {len(files)} datasets")

    # Write to DuckDB
    output_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(output_path))
    con.execute("CREATE TABLE screens AS SELECT * FROM combined")
    con.execute("CREATE INDEX idx_dataset ON screens(Metadata_dataset)")
    con.execute("CREATE INDEX idx_pert_type ON screens(Metadata_pert_type)")
    con.close()

    logger.info(f"Database created at {output_path}")
    return output_path


def _bh_adjusted_critical_value(pvalues: np.ndarray | pd.Series, fdr: float = 0.05) -> float:
    """Calculate Benjamini-Hochberg adjusted critical value for multiple testing.

    Args:
        pvalues: Array of p-values to adjust
        fdr: False discovery rate threshold (default 0.05)

    Returns:
        Adjusted critical value for significance testing
    """
    sorted_pvalues = np.sort(pvalues)
    m = len(pvalues)
    ranks = np.arange(1, m + 1)
    critical_values = (ranks / m) * fdr
    below_threshold = sorted_pvalues <= critical_values
    if np.any(below_threshold):
        adjusted_critical = sorted_pvalues[below_threshold].max()
    else:
        adjusted_critical = np.nan
    return adjusted_critical


def _save_excel_with_sheets(
    filename: Path,
    dataframes: list[pd.DataFrame],
    sheet_names: list[str],
    keep_index: bool = True,
) -> None:
    """Save DataFrames as sheets in Excel file (creates new or overwrites existing).

    Args:
        filename: Path to Excel file
        dataframes: List of DataFrames to save
        sheet_names: List of sheet names (must match length of dataframes)
        keep_index: Whether to keep DataFrame index in Excel output
    """
    if len(dataframes) != len(sheet_names):
        raise ValueError("The number of DataFrames must match the number of sheet names.")

    with pd.ExcelWriter(filename, engine="openpyxl", mode="w") as writer:
        for df, sheet_name in zip(dataframes, sheet_names, strict=True):
            df.to_excel(writer, sheet_name=sheet_name, index=keep_index)


def process_single_virtual_screen_csv(
    dataset: str,
    csv_path: Path,
    output_dir: Path,
    fdr: float = 0.05,
    cell_count_quantile: float | None = None,
    sort_by_col: str = "d_slope",
    p_val_target_col: str = "p_slope_std",
    p_val_orth_col: str = "p_orth_std",
    parquet_output_dir: Path | None = None,
) -> dict[str, int]:
    """Process a single virtual screen CSV file to Excel and optionally Parquet.

    Args:
        dataset: Dataset name (e.g., "CDRP", "lincs", "jump_orf")
        csv_path: Path to input CSV file
        output_dir: Directory to save output Excel file
        fdr: False discovery rate for Benjamini-Hochberg correction
        cell_count_quantile: If set, filter out bottom quantile by cell count (e.g., 0.1 for bottom 10%)
        sort_by_col: Column to sort results by (default: "d_slope")
        p_val_target_col: Column name for target feature p-values
        p_val_orth_col: Column name for orthogonal features p-values
        parquet_output_dir: If set, save unfiltered data as Parquet file to this directory

    Returns:
        Dictionary of filtering statistics with keys:
        - 'raw': Number of perturbations before filtering
        - 'target_sig': After target feature significance filter
        - 'orth_filt': After orthogonal feature filter
        - 'both_filt': After both filters
    """
    # Ensure output directories exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if parquet_output_dir is not None:
        parquet_output_dir = Path(parquet_output_dir)
        parquet_output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    logger.info(f"Processing {dataset}...")

    # Get dataset metadata configuration
    meta_cols = DATASET_INFO[dataset]["meta_cols"]
    pert_col = DATASET_INFO[dataset]["pert_col"]

    # Read CSV
    res_df = pd.read_csv(csv_path)
    logger.info(f"  Loaded {len(res_df):,} rows")

    # Remove rows with null values in sort column
    res_df = res_df[~res_df[sort_by_col].isnull()].reset_index(drop=True)

    # Initialize filtering stats for this dataset
    stats = {"raw": res_df.shape[0]}

    # Apply cell count filter if specified
    if cell_count_quantile is not None:
        threshold = res_df["Count_Cells_avg"].quantile(cell_count_quantile)
        res_df_filtered = res_df[res_df["Count_Cells_avg"] > threshold].reset_index(drop=True)
        logger.info(f"  Cell count filter: {len(res_df_filtered):,} rows (removed bottom {cell_count_quantile*100:.0f}%)")
    else:
        res_df_filtered = res_df.copy()

    # Calculate BH-corrected critical values
    target_bh_critical = np.round(_bh_adjusted_critical_value(res_df_filtered[p_val_target_col], fdr=fdr), 5)
    orth_bh_critical = np.round(_bh_adjusted_critical_value(res_df_filtered[p_val_orth_col], fdr=fdr), 5)

    logger.info(f"  BH critical values: target={target_bh_critical:.5f}, orth={orth_bh_critical:.5f}")

    # Apply target feature significance filter
    res_df_target_sig = res_df_filtered[res_df_filtered[p_val_target_col] < target_bh_critical].reset_index(drop=True)
    stats["target_sig"] = res_df_target_sig.shape[0]

    # Apply orthogonal feature filter
    res_df_orth_filt = res_df_filtered[res_df_filtered[p_val_orth_col] > orth_bh_critical].reset_index(drop=True)
    stats["orth_filt"] = res_df_orth_filt.shape[0]

    # Apply both filters
    res_df_both_filt = res_df_target_sig[res_df_target_sig[p_val_orth_col] > orth_bh_critical].reset_index(drop=True)
    stats["both_filt"] = res_df_both_filt.shape[0]

    # Define columns to save
    cols_to_save = [pert_col, sort_by_col, p_val_target_col, p_val_orth_col, "t_orth", "Count_Cells_avg"] + meta_cols

    # Prepare three versions for saving
    # Sheet 1: All data (sorted)
    df_all = res_df.sort_values(by=sort_by_col)[cols_to_save]
    df_all = df_all.loc[:, ~df_all.columns.duplicated(keep="last")]

    # Sheet 2: Orthogonal filter only (sorted)
    df_orth = res_df_orth_filt.sort_values(by=sort_by_col)[cols_to_save]
    df_orth = df_orth.loc[:, ~df_orth.columns.duplicated(keep="last")]

    # Sheet 3: Both filters (sorted)
    df_both = res_df_both_filt.sort_values(by=sort_by_col)[cols_to_save]
    df_both = df_both.loc[:, ~df_both.columns.duplicated(keep="last")]

    # Save to Excel
    output_path = output_dir / f"{dataset}_screen_results.xlsx"
    _save_excel_with_sheets(
        output_path,
        [df_all, df_orth, df_both],
        [dataset, f"{dataset}_orthfilt", f"{dataset}_bothfilt"],
        keep_index=False,
    )

    logger.info(f"  Saved Excel to: {output_path}")

    # Save unfiltered data to Parquet if requested
    if parquet_output_dir is not None:
        parquet_path = parquet_output_dir / f"{dataset}_unfiltered.parquet"
        df_all.to_parquet(parquet_path, index=False)
        logger.info(f"  Saved Parquet to: {parquet_path}")

    logger.info(f"  Filtering: raw={stats['raw']} -> target_sig={stats['target_sig']} -> orth_filt={stats['orth_filt']} -> both_filt={stats['both_filt']}")

    return stats


def validate_databases(baseline_path: Path, new_path: Path) -> bool:
    """Compare two DuckDB databases for identical data.

    Validates that two DuckDB databases contain identical data by comparing
    all rows and columns, accounting for NaN values properly.

    Args:
        baseline_path: Path to baseline DuckDB file
        new_path: Path to new DuckDB file to validate

    Returns:
        True if databases are identical, False otherwise
    """
    baseline_path = Path(baseline_path)
    new_path = Path(new_path)

    if not baseline_path.exists():
        logger.error(f"Baseline file not found: {baseline_path}")
        return False

    if not new_path.exists():
        logger.error(f"New file not found: {new_path}")
        return False

    logger.info(f"Comparing databases:")
    logger.info(f"  Baseline: {baseline_path}")
    logger.info(f"  New:      {new_path}")

    # Connect to both databases
    baseline = duckdb.connect(str(baseline_path), read_only=True)
    new = duckdb.connect(str(new_path), read_only=True)

    # Load data sorted for consistent comparison
    df_baseline = baseline.execute("SELECT * FROM screens ORDER BY Metadata_dataset, Metadata_pert_id").df()
    df_new = new.execute("SELECT * FROM screens ORDER BY Metadata_dataset, Metadata_pert_id").df()

    baseline.close()
    new.close()

    logger.info(f"Baseline rows: {len(df_baseline):,}")
    logger.info(f"New rows: {len(df_new):,}")

    # Compare DataFrames accounting for NaN values
    identical = True

    if df_baseline.shape != df_new.shape:
        logger.error(f"Shape differs: baseline={df_baseline.shape}, new={df_new.shape}")
        identical = False
    elif list(df_baseline.columns) != list(df_new.columns):
        logger.error(f"Columns differ")
        logger.error(f"  Baseline: {list(df_baseline.columns)}")
        logger.error(f"  New: {list(df_new.columns)}")
        identical = False
    else:
        # Check each column, accounting for NaN values
        for col in df_baseline.columns:
            if df_baseline[col].dtype == "object":
                # For object columns, compare with NaN handling
                mask = ~(df_baseline[col].isna() & df_new[col].isna())
                diff = (df_baseline.loc[mask, col] != df_new.loc[mask, col]).any()
            else:
                # For numeric columns, use np.allclose
                diff = not np.allclose(df_baseline[col], df_new[col], equal_nan=True)

            if diff:
                logger.error(f"Column {col} differs")
                identical = False
                break

    if identical:
        logger.info("✓ SUCCESS: Databases are identical")
    else:
        logger.error("✗ FAILURE: Databases differ")

    return identical
