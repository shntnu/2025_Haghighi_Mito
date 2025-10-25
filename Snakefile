"""Snakemake pipeline for mitochondrial morphology screen analysis."""

# Default rule - runs all pipeline steps
rule all:
    input:
        "data/processed/screen_results.duckdb"


# Convert virtual screen CSV results to filtered Excel tables
rule process_virtual_screen_results:
    input:
        csv_dir="data/external/mito_project/workspace/results/virtual_screen",
        csvs=[
            "data/external/mito_project/workspace/results/virtual_screen/CDRP_results_pattern_aug_070624.csv",
            "data/external/mito_project/workspace/results/virtual_screen/jump_compound_results_pattern_aug_070624.csv",
            "data/external/mito_project/workspace/results/virtual_screen/jump_crispr_results_pattern_aug_070624.csv",
            "data/external/mito_project/workspace/results/virtual_screen/jump_orf_results_pattern_aug_070624.csv",
            "data/external/mito_project/workspace/results/virtual_screen/lincs_results_pattern_aug_070624.csv",
            "data/external/mito_project/workspace/results/virtual_screen/taorf_results_pattern_aug_070624.csv",
        ]
    output:
        "data/processed/tables/generated_from_s3_baseline/CDRP_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_compound_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_crispr_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_orf_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/lincs_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/taorf_screen_results.xlsx"
    params:
        output_dir="data/processed/tables/generated_from_s3_baseline"
    shell:
        "pixi run haghighi-mito process-csvs --output-dir {params.output_dir}"


# Create unified DuckDB database from individual Excel screen results
rule create_database:
    input:
        "data/processed/tables/generated_from_s3_baseline/CDRP_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_compound_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_crispr_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/jump_orf_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/lincs_screen_results.xlsx",
        "data/processed/tables/generated_from_s3_baseline/taorf_screen_results.xlsx"
    output:
        "data/processed/screen_results.duckdb"
    params:
        tables_dir="data/processed/tables/generated_from_s3_baseline",
        output_path="data/processed/screen_results.duckdb"
    shell:
        "pixi run haghighi-mito create-database --output-path {params.output_path} --tables-dir {params.tables_dir} --overwrite"
