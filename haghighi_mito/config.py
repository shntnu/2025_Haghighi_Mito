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

# External curated reference tables (manually curated, not pipeline outputs)
EXTERNAL_TABLES_DIR = EXTERNAL_DATA_DIR / "tables" / "curated_2025-10-25"
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
FIBROBLAST_SINGLECELL_PARQUET = INTERIM_DATA_DIR / "fibroblast_singlecell.parquet"
PATIENT_LABELS_PATH = EXTERNAL_TABLES_DIR / "patient_labels_updatedSept302025.csv"
PATIENT_LABELS_FULL_PATH = EXTERNAL_TABLES_DIR / "patient_labels_updatedSept302025_full.csv"
PATIENT_FIGURES_DIR = PROCESSED_DATA_DIR / "figures" / "patient_phenotype"
SUPPLEMENTAL_FIGURES_DIR = PROCESSED_DATA_DIR / "figures" / "supplemental"

# Patient category display constants
PATIENT_ORDER = ["Control", "psychosis", "BP", "SZ", "SZA", "MDD"]
PATIENT_PALETTE = ["lightgray", "#a484ac", "firebrick", "lightcoral", "pink", "lightsteelblue"]

# ---------------------------------------------------------------------------
# Imaging calibration: pixel size of the patient fibroblast images
# ---------------------------------------------------------------------------
# Used by notebooks/3.0-mh-slope-analysis.py to convert CellProfiler AreaShape
# features (reported in pixels) into microns for Supp Fig 3 panels A–E.
#
# !!! IMPORTANT: THIS VALUE IS *INFERRED*, NOT MEASURED. !!!
#
#   - No calibration slide / stage micrometer image exists in the repo or in
#     the published Zenodo data (doi:10.5281/zenodo.15390513), so the pixel
#     size cannot be ground-truth verified.
#   - The value below is derived from three indirect pieces of evidence
#     (protocol text, TIFF pixel dimensions, and a sanity check against
#     expected fibroblast cell sizes). The camera model is not recorded in
#     the acquisition protocol — it is deduced from the image dimensions
#     matching one specific Zeiss camera's native sensor resolution.
#   - Until confirmed by the person who collected the data (Erin Weisbart),
#     treat raw-micron figures as provisional. Changing this single constant
#     and re-running `snakemake generate_slope_figures` will regenerate all
#     raw-micron outputs. Welch t-test p-values and Pearson correlations are
#     affine-invariant, so statistical conclusions do not depend on the
#     exact value.
#
# ===========================================================================
# Evidence chain for PIXEL_SIZE_UM = 0.161 µm/px (inferred, matches upstream)
# ===========================================================================
#
# 1) PROTOCOL (protocols/McleanCollectionFibroblastGrowthProtocol.md:122)
#    Quote: "25 fields from each of the two coverslips (50 fields per sample)
#    were imaged with a 40X objective in three independent channels on a Zeiss
#    Axiovert Observer 2.1 with a Colibri 2 LED illumination system. Image
#    acquisition was done with Zeiss Zen 2.0 software."
#    → Objective magnification = 40× (confirmed).
#    → Acquisition software = Zeiss ZEN (confirmed in TIFF tag 305 below).
#    → Camera model is NOT stated — has to be inferred from the image files.
#
# 2) TIFF DIMENSIONS (all 24 source images identical)
#    Inspecting data/processed/figures/patient_phenotype/figure2_source_images/
#    all 24 TIFFs are:
#        shape      = (1040, 1388)         # height × width, 1388 × 1040 px
#        dtype      = uint16                # 12–14 bit monochrome, stored 16-bit
#        Software   = "ZEN 2012 (blue edition)"   (tag 305)
#        XResolution = 96000/1000 = 96 dpi  (tag 282, ResolutionUnit=INCH)
#    → 1388 × 1040 is the native sensor resolution of exactly one Zeiss camera
#      sold with the Axiovert Observer 2.1: the AxioCam MRm (monochrome).
#      No other AxioCam model has this exact pixel grid. See Zeiss datasheet
#      "AxioCam MRm / MRc", 1388 × 1040 px, 6.45 µm × 6.45 µm pixel pitch,
#      8.95 mm × 6.71 mm sensor area (2/3" format, Sony ICX285 CCD).
#    → The 96 dpi XResolution is a Windows print-resolution placeholder that
#      ZEN writes on TIFF export; it is NOT physical calibration. Real
#      calibration ("ScalingX/Y" in meters/pixel) lives only in the native
#      .czi files, which are stripped on TIFF export.
#    → Matches the date evidence: TIFF XMP "CreateDate=2014-11-20"; AxioCam MRm
#      was the standard monochrome fluorescence camera sold with the Observer
#      platform in the 2012–2015 timeframe.
#
# 3) ARITHMETIC
#        pixel_size = sensor_pitch / objective_magnification
#                   = 6.45 µm / 40
#                   = 0.16125 µm/px  →  area factor ≈ 0.026 µm²/px²
#
#    Upstream (Erin's plot_singlecell_cellsize.ipynb, cell 14) uses the
#    3-sig-fig value ``microns_per_pixel = 0.161``. We adopt the same value
#    so figures are numerically identical to the published Supp Fig 3 panels.
#    The 0.155% difference vs the camera-exact 0.16125 is cosmetic and below
#    plot resolution.
#
# 4) SANITY CHECK against the data (SQL-aggregated means per subject)
#    After converting with 0.161 µm/px, whole-cell means land at:
#        Area        ≈ 2,800 µm²    (per-patient means, post-QC 39017 cells)
#        MajorAxis   ≈ 80 µm
#        MinorAxis   ≈ 40 µm
#        Perimeter   ≈ 300 µm
#    These are in the expected range for cultured human dermal fibroblasts
#    (typical area 1,000–3,000 µm²; major-axis 50–150 µm). A wrong pixel size
#    (e.g., a 20× value of 0.3225 µm/px) would put areas at ~11,300 µm² — far
#    too large. A 63× value (0.102 µm/px) would put areas at ~1,130 µm² — on
#    the low end. The 40× assumption produces physically plausible cells.
#
# ===========================================================================
# Caveats
# ===========================================================================
# - No binning evidence: the protocol does not mention binning, and 1388 × 1040
#   is the un-binned native resolution, so 1×1 binning is assumed.
# - No calibration slide (e.g., a stage micrometer image) is shipped with the
#   repo, so this cannot be ground-truth verified from the data alone.
# - If Erin confirms the acquisition used a different camera or binning,
#   updating this single constant and re-running the figure rule is enough —
#   relative values and p-values are unchanged (affine-invariant).
PIXEL_SIZE_UM = 0.161

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
