# Mitochondrial Morphology Screen Pipeline
# Lightweight task runner for Snakemake pipeline
#
# See Snakefile docstring for comprehensive pipeline documentation.
# Three methods available:
#   Method 0: Baseline (validated, July 2024)
#   Method 1: Notebook (complete but messy, 1433 lines)
#   Method 2: Clean module (incomplete, 448 lines)

# Show available commands
default:
    @just --list

# ============================================================================
# METHOD 0: BASELINE PIPELINE (Recommended for Production)
# ============================================================================

# Download pre-computed baseline results from S3 (65 MB, validated July 2024)
download-baseline:
    pixi run snakemake download_all_baseline --cores 4 --printshellcmds

# Run baseline pipeline (process S3 CSVs → Excel + DuckDB, ~5 min)
run-baseline:
    pixi run snakemake all_baseline --cores 4 --printshellcmds

# ============================================================================
# SHARED: Download Raw Data for Regeneration (Methods 1 & 2)
# ============================================================================

# Download raw data for regenerating results (2.7 GB: per-site profiles + metadata)
download-raw:
    pixi run snakemake download_screening_data --cores 4 --printshellcmds

# ============================================================================
# METHOD 1: REGENERATED - Notebook (Complete Pipeline)
# ============================================================================

# [Method 1] Run full notebook pipeline for all datasets (CSV → Excel → DuckDB, ~10 min/dataset)
run-notebook:
    pixi run snakemake all_notebook --cores 4 --printshellcmds

# [Method 1] Run notebook for a specific dataset (CSV only, use run-notebook for full pipeline)
run-notebook-for DATASET:
    pixi run snakemake data/external/mito_project/workspace/results/virtual_screen_notebook/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# ============================================================================
# METHOD 2: REGENERATED - Clean Module (Incomplete - Stops at CSV)
# ============================================================================

# TODO: Add when Method 2 Excel/DuckDB processing is complete
# [Method 2] Run full module pipeline for all datasets (CSV → Excel → DuckDB)
# run-module:
#     pixi run snakemake all_module --cores 4 --printshellcmds

# [Method 2] Run clean module for a specific dataset (CSV + comparison only, ⚠️ no Excel/DuckDB)
run-module-for DATASET:
    pixi run snakemake data/processed/virtual_screen_module/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# [Method 2] Generate baseline comparison plots for a specific dataset
plot-comparison-for DATASET:
    pixi run snakemake data/processed/figures/t_target_pattern_analysis/{{DATASET}}_baseline_vs_regenerated.png --cores 1 --printshellcmds

# [Method 2] Run clean module for all datasets (CSV + comparison only, ⚠️ no Excel/DuckDB)
run-all-modules:
    pixi run snakemake run_all_virtual_screen_modules --cores 4 --printshellcmds

# [Method 2] Generate all baseline comparison plots
plot-all-comparisons:
    pixi run snakemake plot_all_baseline_comparisons --cores 4 --printshellcmds

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

# Clean baseline outputs (Method 0)
clean-baseline:
    rm -f data/processed/screen_results_baseline.duckdb
    rm -rf data/interim/parquet_baseline/
    rm -rf data/processed/tables/generated_from_s3_baseline/

# Clean notebook outputs (Method 1)
clean-notebook:
    rm -f data/processed/screen_results_notebook.duckdb
    rm -rf data/interim/parquet_notebook/
    rm -rf data/processed/tables/generated_from_notebook/

# Clean module outputs (Method 2)
clean-module:
    rm -f data/processed/screen_results_module.duckdb
    rm -rf data/processed/virtual_screen_module/
    rm -rf data/processed/figures/t_target_pattern_analysis/
    rm -rf data/interim/parquet_module/
    rm -rf data/processed/tables/generated_from_module/

# Clean everything
clean-all: clean-baseline clean-notebook clean-module
