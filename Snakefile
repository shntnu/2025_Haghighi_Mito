"""Snakemake pipeline for mitochondrial morphology screen analysis."""

# Default rule - runs all pipeline steps
rule all:
    input:
        "data/processed/screen_results.duckdb"


# Create unified DuckDB database from individual Excel screen results
rule create_database:
    input:
        "data/processed/tables/CDRP_screen_results.xlsx",
        "data/processed/tables/jump_compound_screen_results.xlsx",
        "data/processed/tables/jump_crispr_screen_results.xlsx",
        "data/processed/tables/jump_orf_screen_results.xlsx",
        "data/processed/tables/lincs_screen_results.xlsx",
        "data/processed/tables/taorf_screen_results.xlsx"
    output:
        "data/processed/screen_results.duckdb"
    shell:
        "uv run python -c 'from haghighi_mito.data import create_screen_database; create_screen_database(overwrite=True)'"
