"""Project configuration and path management."""

from pathlib import Path

from dotenv import load_dotenv

# Load environment variables from .env file if it exists
load_dotenv()

# Paths
PROJ_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"

# Processed data subdirectories
PROCESSED_TABLES_DIR = PROCESSED_DATA_DIR / "tables" / "curated_2025-10-25"
PROCESSED_FIGURES_DIR = PROCESSED_DATA_DIR / "figures"
