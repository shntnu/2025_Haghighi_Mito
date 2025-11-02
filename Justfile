# Mitochondrial Morphology Screen Pipeline
# Lightweight task runner for Snakemake pipeline
#
# See Snakefile docstring for comprehensive pipeline documentation.
# Three methods available:
#   Method 0: Baseline (validated, July 2024)
#   Method 1: Notebook (original exploratory implementation)
#   Method 2: Module (clean refactored implementation)

# Show available commands
default:
    @just --list --unsorted

# ============================================================================
# METHOD 0: BASELINE PIPELINE (Recommended for Production)
# ============================================================================

# (metrics CSV → Excel/Parquet → DuckDB)
generate-baseline-all:
    pixi run snakemake all_baseline --cores 4 --printshellcmds

# ============================================================================
# METHOD 1: REGENERATED - Notebook (Original Implementation)
# ============================================================================

# [Method 1] (profiles → metrics CSV → Excel/Parquet → DuckDB)
generate-notebook-all:
    pixi run snakemake all_notebook --cores 4 --printshellcmds

# [Method 1] (profiles → metrics CSV)
generate-notebook-csvs:
    pixi run snakemake all_notebook_csvs --cores 4 --printshellcmds

# [Method 1] Single dataset metrics CSV
generate-notebook-csv-for DATASET:
    pixi run snakemake data/processed/virtual_screen_notebook/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# ============================================================================
# METHOD 2: REGENERATED - Module (Refactored Implementation)
# ============================================================================

# [Method 2] (profiles → metrics CSV → Excel/Parquet → DuckDB)
generate-module-all:
    pixi run snakemake all_module --cores 4 --printshellcmds

# [Method 2] (profiles → metrics CSV)
generate-module-csvs:
    pixi run snakemake all_module_csvs --cores 4 --printshellcmds

# [Method 2] Single dataset metrics CSV
generate-module-csv-for DATASET:
    pixi run snakemake data/processed/virtual_screen_module/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# ============================================================================
# DIAGNOSTICS & DEBUGGING
# ============================================================================

# Compare module results with baseline (all datasets, fast)
diagnose-all:
    pixi run snakemake all_module_diagnose --cores 4 --printshellcmds

# Compare module results with baseline (single dataset, fast)
diagnose-for DATASET:
    pixi run snakemake data/processed/virtual_screen_module/{{DATASET}}_comparison_metrics.png --cores 1 --printshellcmds

# Recalculate slopes from raw data (slow, educational alternative to diagnose-for)
reproduce-for DATASET="taorf":
    @echo "Downloading baseline, per-site data, and metadata for {{DATASET}} if needed..."
    pixi run snakemake \
        data/external/mito_project/workspace/results/virtual_screen_baseline/{{DATASET}}_results_pattern_aug_070624.csv \
        data/external/mito_project/workspace/per_site_aggregated_profiles_newpattern_2/{{DATASET}}/.download_complete \
        data/external/mito_project/workspace/metadata/CDRP_meta.csv \
        data/external/mito_project/workspace/metadata/JUMP-ORF/ORF_list.tsv \
        data/external/mito_project/workspace/metadata/JUMP/compound.csv.gz \
        data/external/mito_project/workspace/metadata/JUMP/crispr.csv.gz \
        data/external/mito_project/workspace/metadata/JUMP/orf.csv.gz \
        data/external/mito_project/workspace/metadata/JUMP/plate.csv.gz \
        data/external/mito_project/workspace/metadata/JUMP/well.csv.gz \
        data/external/mito_project/workspace/metadata/LINCS_meta.csv \
        data/external/mito_project/workspace/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz \
        data/external/mito_project/workspace/metadata/lincs/DrugRepurposing_Metadata.csv \
        --cores 1 --printshellcmds
    @echo "Running slope discrepancy analysis..."
    pixi run python scripts/reproduce_slope_discrepancy.py {{DATASET}}

# ============================================================================
# UTILITIES
# ============================================================================

# Generate DAG visualizations
viz:
    ./scripts/generate-dag.sh

# Preview what will run (dry run)
dry:
    pixi run snakemake -n -p

# Show pipeline status
status:
    pixi run snakemake --summary

# Compare databases (e.g., baseline vs notebook)
validate-databases DB1 DB2:
    pixi run haghighi-mito validate-databases \
        --baseline data/processed/screen_results_{{DB1}}.duckdb \
        --new data/processed/screen_results_{{DB2}}.duckdb

