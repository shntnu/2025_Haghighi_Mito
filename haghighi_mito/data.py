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

    # Screen result files
    files = {
        "CDRP": PROCESSED_TABLES_DIR / "CDRP_screen_results.xlsx",
        "JUMP_Compound": PROCESSED_TABLES_DIR / "jump_compound_screen_results.xlsx",
        "JUMP_CRISPR": PROCESSED_TABLES_DIR / "jump_crispr_screen_results.xlsx",
        "JUMP_ORF": PROCESSED_TABLES_DIR / "jump_orf_screen_results.xlsx",
        "LINCS": PROCESSED_TABLES_DIR / "lincs_screen_results.xlsx",
        "TA_ORF": PROCESSED_TABLES_DIR / "taorf_screen_results.xlsx",
    }

    # Load and combine all datasets
    dfs = []
    for dataset_name, file_path in files.items():
        logger.info(f"Loading {dataset_name} from {file_path.name}")
        df = pd.read_excel(file_path)
        df.insert(0, "dataset_source", dataset_name)
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined {len(combined):,} rows from {len(files)} datasets")

    # Write to DuckDB
    output_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(output_path))
    con.execute("CREATE TABLE screens AS SELECT * FROM combined")
    con.execute("CREATE INDEX idx_dataset ON screens(dataset_source)")
    con.close()

    logger.info(f"Database created at {output_path}")
    return output_path
