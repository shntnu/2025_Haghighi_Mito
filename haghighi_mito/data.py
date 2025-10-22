"""Data loading and database utilities."""

from pathlib import Path

import duckdb
import pandas as pd
from loguru import logger

from haghighi_mito.config import PROCESSED_DATA_DIR, PROCESSED_TABLES_DIR


def create_screen_database(
    output_path: Path = PROCESSED_DATA_DIR / "screen_results.duckdb",
    overwrite: bool = False,
) -> Path:
    """Combine all screen result Excel files into a single DuckDB database.

    Only relevant columns are kept for each dataset to reduce database size.

    Args:
        output_path: Path to output DuckDB database file
        overwrite: If True, recreate database even if it exists

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

    # Screen result files
    files = {
        "CDRP": PROCESSED_TABLES_DIR / "CDRP_screen_results.xlsx",
        "JUMP_Compound": PROCESSED_TABLES_DIR / "jump_compound_screen_results.xlsx",
        "JUMP_CRISPR": PROCESSED_TABLES_DIR / "jump_crispr_screen_results.xlsx",
        "JUMP_ORF": PROCESSED_TABLES_DIR / "jump_orf_screen_results.xlsx",
        "LINCS": PROCESSED_TABLES_DIR / "lincs_screen_results.xlsx",
        "TA_ORF": PROCESSED_TABLES_DIR / "taorf_screen_results.xlsx",
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
        logger.info(f"Loading {dataset_name} from {file_path.name}")
        df = pd.read_excel(file_path)

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
