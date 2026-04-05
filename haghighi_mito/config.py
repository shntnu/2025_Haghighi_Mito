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

# Mito project paths (from S3 download)
MITO_PROJECT_DIR = EXTERNAL_DATA_DIR / "mito_project"
MITO_WORKSPACE_DIR = MITO_PROJECT_DIR / "workspace"
MITO_VIRTUAL_SCREEN_DIR = MITO_WORKSPACE_DIR / "results" / "virtual_screen"
MITO_ORTH_FEATURES_DIR = MITO_WORKSPACE_DIR / "results" / "target_pattern_orth_features_lists"

# Patient fibroblast analysis paths
FIBROBLAST_DATA_DIR = MITO_WORKSPACE_DIR / "singleCellData"
SQLITE_DATA_DIR = MITO_WORKSPACE_DIR / "backend" / "Mito_Morphology_input"
AGGREGATED_PROFILES_PATH = FIBROBLAST_DATA_DIR / "aggregated_profiles_fibroblast.csv"
PATIENT_LABELS_PATH = PROCESSED_TABLES_DIR / "patient_labels_updatedSept302025.csv"
PATIENT_LABELS_FULL_PATH = PROCESSED_TABLES_DIR / "patient_labels_updatedSept302025_full.csv"
PATIENT_FIGURES_DIR = PROCESSED_DATA_DIR / "figures" / "patient_phenotype"
SUPPLEMENTAL_FIGURES_DIR = PROCESSED_DATA_DIR / "figures" / "supplemental"

# Patient category display constants
PATIENT_ORDER = ["Control", "psychosis", "BP", "SZ", "SZA", "MDD"]
PATIENT_PALETTE = ["lightgray", "#a484ac", "firebrick", "lightcoral", "pink", "lightsteelblue"]

# Dataset metadata configuration for virtual screen analysis
DATASET_INFO = {
    "taorf": {
        "meta_cols": [
            "Metadata_gene_name",
            "Metadata_pert_name",
            "Metadata_broad_sample",
            "Metadata_moa",
        ],
        "pert_col": "Metadata_broad_sample",
    },
    "CDRP": {
        "meta_cols": [
            "Metadata_broad_sample",
            "Metadata_mmoles_per_liter2",
            "Metadata_pert_id",
            "Metadata_Sample_Dose",
            "Metadata_moa",
        ],
        "pert_col": "Metadata_Sample_Dose",
    },
    "lincs": {
        "meta_cols": [
            "Metadata_broad_sample",
            "Metadata_dose_recode",
            "Metadata_pert_id",
            "Metadata_pert_mfc_id",
            "Metadata_InChIKey14",
            "Metadata_pert_type",
            "Metadata_moa",
            "Metadata_target",
            "Metadata_pert_id_dose",
            "Metadata_pert_name",
        ],
        "pert_col": "Metadata_pert_id_dose",
    },
    "jump_orf": {
        "meta_cols": ["Metadata_Symbol", "Metadata_broad_sample", "Metadata_JCP2022"],
        "pert_col": "Metadata_JCP2022",
    },
    "jump_crispr": {
        "meta_cols": ["Metadata_NCBI_Gene_ID", "Metadata_Symbol", "Metadata_JCP2022"],
        "pert_col": "Metadata_JCP2022",
    },
    "jump_compound": {
        "meta_cols": ["Metadata_InChIKey", "Metadata_InChI", "Metadata_JCP2022"],
        "pert_col": "Metadata_JCP2022",
    },
}
