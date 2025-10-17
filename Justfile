# Mitochondrial Morphology Screen Pipeline
# Usage: just <command>

# Compute resources
CORES := env_var_or_default("CORES", "4")

# Snakemake configuration
SNAKEMAKE := "uv run snakemake"

# Show available commands
default:
    @just --list

# Run the full pipeline (create database from Excel files)
run:
    @echo "Running pipeline with {{CORES}} cores..."
    @{{SNAKEMAKE}} --cores {{CORES}} --printshellcmds

# Preview what will run without executing
dry:
    @echo "Performing dry run..."
    @{{SNAKEMAKE}} --cores {{CORES}} -n -p

# Delete generated database
clean:
    @echo "Cleaning generated database..."
    @rm -f data/processed/screen_results.duckdb

# Show pipeline status
status:
    @{{SNAKEMAKE}} --summary

# Show current configuration
config:
    @echo "Current Configuration:"
    @echo "  CORES: {{CORES}}"
    @echo "  Database: data/processed/screen_results.duckdb"
