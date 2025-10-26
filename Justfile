# Mitochondrial Morphology Screen Pipeline
# Lightweight task runner for Snakemake pipeline

# Show available commands
default:
    @just --list

# Run baseline pipeline (process S3 results → Excel + DuckDB)
run-baseline:
    pixi run snakemake all_baseline --cores 4 --printshellcmds

# Run regenerated pipeline (process locally-regenerated results → Excel + DuckDB)
run-regenerated:
    pixi run snakemake all_regenerated --cores 4 --printshellcmds

# Download pre-computed baseline results from S3 (6 CSVs, for processing pipeline)
download-baseline:
    pixi run snakemake download_all_baseline --cores 4 --printshellcmds

# Download raw data for regenerating results (per-site profiles + metadata, for notebook 2.0)
download-raw:
    pixi run snakemake download_screening_data --cores 4 --printshellcmds

# Run notebook for a specific dataset (e.g., just run-notebook-for taorf)
run-notebook-for DATASET:
    pixi run snakemake data/external/mito_project/workspace/results/virtual_screen_regenerated/{{DATASET}}_results_pattern_aug_070624.csv --cores 1 --printshellcmds

# Run virtual screen analysis for a specific dataset (e.g., just run-virtual-screen-for taorf)
run-virtual-screen-for DATASET:
    pixi run haghighi-mito virtual-screen --dataset {{DATASET}} --compare-baseline

# Create baseline comparison plots for a specific dataset (e.g., just plot-baseline-comparison-for taorf)
plot-baseline-comparison-for DATASET:
    pixi run haghighi-mito plot-baseline-comparison --dataset {{DATASET}}

# Generate DAG visualizations (simplified rules + full jobs for both pipelines)
viz:
    #!/usr/bin/env bash
    mkdir -p docs/pipeline

    # Generate simplified DAGs (rules only)
    pixi run snakemake all_baseline --rulegraph 2>&1 | tail -n +2 > /tmp/rg_baseline.dot
    pixi run snakemake all_regenerated --rulegraph 2>&1 | tail -n +2 > /tmp/rg_regenerated.dot
    dot -Tpng /tmp/rg_baseline.dot -o docs/pipeline/dag_baseline_rules.png
    dot -Tpng /tmp/rg_regenerated.dot -o docs/pipeline/dag_regenerated_rules.png

    # Generate full DAGs (all job instances)
    pixi run snakemake all_baseline --dag 2>&1 | tail -n +2 > /tmp/dag_baseline.dot
    pixi run snakemake all_regenerated --dag 2>&1 | tail -n +2 > /tmp/dag_regenerated.dot
    dot -Tpng /tmp/dag_baseline.dot -o docs/pipeline/dag_baseline_jobs.png
    dot -Tpng /tmp/dag_regenerated.dot -o docs/pipeline/dag_regenerated_jobs.png

    # Cleanup temp files
    rm /tmp/rg_baseline.dot /tmp/rg_regenerated.dot /tmp/dag_baseline.dot /tmp/dag_regenerated.dot

    echo "Generated in docs/pipeline/:"
    echo "  - dag_baseline_rules.png (simplified)"
    echo "  - dag_regenerated_rules.png (simplified)"
    echo "  - dag_baseline_jobs.png (full)"
    echo "  - dag_regenerated_jobs.png (full)"

# Preview what will run (dry run)
dry:
    pixi run snakemake -n -p

# Show pipeline status
status:
    pixi run snakemake --summary

# Clean baseline outputs
clean-baseline:
    rm -f data/processed/screen_results_baseline.duckdb
    rm -rf data/interim/parquet_baseline/
    rm -rf data/processed/tables/generated_from_s3_baseline/

# Clean regenerated outputs
clean-regenerated:
    rm -f data/processed/screen_results_regenerated.duckdb
    rm -rf data/interim/parquet_regenerated/
    rm -rf data/processed/tables/generated_from_local/

# Clean everything
clean-all: clean-baseline clean-regenerated
