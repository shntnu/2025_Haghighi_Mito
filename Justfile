# Mitochondrial Morphology Screen Pipeline
# Lightweight task runner for Snakemake pipeline

# Show available commands
default:
    @just --list

# Run baseline pipeline (default)
run:
    pixi run snakemake --cores 4 --printshellcmds

# Run regenerated pipeline
run-regen:
    pixi run snakemake all_regenerated --cores 4 --printshellcmds

# Download baseline CSVs from S3
download:
    pixi run snakemake download_all_baseline --cores 4 --printshellcmds

# Preview what will run (dry run)
dry:
    pixi run snakemake -n -p

# Show pipeline status
status:
    pixi run snakemake --summary

# Clean baseline outputs
clean:
    rm -f data/processed/screen_results_baseline.duckdb
    rm -rf data/interim/parquet_baseline/
    rm -rf data/processed/tables/generated_from_s3_baseline/

# Clean regenerated outputs
clean-regen:
    rm -f data/processed/screen_results_regenerated.duckdb
    rm -rf data/interim/parquet_regenerated/
    rm -rf data/processed/tables/generated_from_local/

# Clean everything
clean-all: clean clean-regen
