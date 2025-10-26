# Mitochondrial Morphology Screen Pipeline
# Lightweight task runner for Snakemake pipeline
#
# See Snakefile docstring for comprehensive pipeline documentation.
# Three methods available:
#   Method 0: Baseline (validated, July 2024)
#   Method 1: Notebook (complete but messy, 1433 lines)
#   Method 2: Clean module (complete, 448 lines)

# Show available commands
default:
    @just --list

# ============================================================================
# METHOD 0: BASELINE PIPELINE (Recommended for Production)
# ============================================================================

# Generate complete baseline pipeline (downloads from S3 if needed, then processes → Excel + DuckDB, ~5 min)
generate-baseline-all:
    pixi run snakemake all_baseline --cores 4 --printshellcmds

# ============================================================================
# METHOD 1: REGENERATED - Notebook (Complete Pipeline)
# ============================================================================

# [Method 1] Generate complete notebook pipeline for all datasets (downloads raw data if needed, then CSV → Excel → DuckDB, ~10 min/dataset)
generate-notebook-all:
    pixi run snakemake all_notebook --cores 4 --printshellcmds

# [Method 1] Generate results CSV for a specific dataset (CSV only, use generate-notebook-all for full pipeline)
generate-notebook-csv-for DATASET:
    pixi run snakemake data/external/mito_project/workspace/results/virtual_screen_notebook/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# [Method 1] Generate results CSVs for all datasets (notebook analysis only, ⚠️ no Excel/DuckDB)
generate-notebook-csvs:
    pixi run snakemake all_notebook_csvs --cores 4 --printshellcmds

# ============================================================================
# METHOD 2: REGENERATED - Clean Module (Complete Pipeline)
# ============================================================================

# [Method 2] Generate complete module pipeline for all datasets (downloads raw data if needed, then CSV → Excel → DuckDB)
generate-module-all:
    pixi run snakemake all_module --cores 4 --printshellcmds

# [Method 2] Generate results CSV for a specific dataset (CSV generation only, ~5-10 min)
generate-module-csv-for DATASET:
    pixi run snakemake data/processed/virtual_screen_module/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# [Method 2] Run diagnostics: compare with baseline + generate plots (FAST - ~1 sec)
diagnose-for DATASET:
    pixi run snakemake data/processed/figures/diagnostics/{{DATASET}}_comparison_metrics.png --cores 1 --printshellcmds

# [Method 2] Generate results CSVs for all datasets (virtual screen analysis only, ⚠️ no Excel/DuckDB)
generate-module-csvs:
    pixi run snakemake all_module_csvs --cores 4 --printshellcmds

# [Method 2] Run diagnostics for all datasets (comparison CSV + plots)
diagnose-all:
    pixi run snakemake all_module_diagnostics --cores 4 --printshellcmds

# ============================================================================
# UTILITIES
# ============================================================================

# Generate DAG visualizations (simplified rules + full jobs for both pipelines)
viz:
    ./scripts/generate-dag.sh

# Preview what will run (dry run)
dry:
    pixi run snakemake -n -p

# Show pipeline status
status:
    pixi run snakemake --summary

# Compare two DuckDB databases for validation (e.g., 'baseline' vs 'notebook')
validate-databases DB1 DB2:
    pixi run haghighi-mito validate-databases \
        --baseline data/processed/screen_results_{{DB1}}.duckdb \
        --new data/processed/screen_results_{{DB2}}.duckdb

