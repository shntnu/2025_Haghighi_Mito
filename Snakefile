"""Snakemake pipeline for mitochondrial morphology screen analysis."""

# Dataset configuration
DATASETS = ["CDRP", "jump_compound", "jump_crispr", "jump_orf", "lincs", "taorf"]

# Default rule - runs all pipeline steps
rule all:
    input:
        "data/processed/screen_results.duckdb"


# Process a single virtual screen CSV file to Excel and Parquet
rule process_single_csv:
    input:
        csv="data/external/mito_project/workspace/results/virtual_screen/{dataset}_results_pattern_aug_070624.csv"
    output:
        excel="data/processed/tables/generated_from_s3_baseline/{dataset}_screen_results.xlsx",
        parquet="data/interim/parquet/{dataset}_unfiltered.parquet"
    params:
        excel_dir="data/processed/tables/generated_from_s3_baseline",
        parquet_dir="data/interim/parquet"
    shell:
        """
        pixi run haghighi-mito process-csv-single \
            --dataset {wildcards.dataset} \
            --csv-path {input.csv} \
            --output-dir {params.excel_dir} \
            --parquet-output-dir {params.parquet_dir}
        """


# Create unified DuckDB database from Parquet files
rule create_database:
    input:
        expand("data/interim/parquet/{dataset}_unfiltered.parquet", dataset=DATASETS)
    output:
        "data/processed/screen_results.duckdb"
    params:
        parquet_dir="data/interim/parquet",
        output_path="data/processed/screen_results.duckdb"
    shell:
        "pixi run haghighi-mito create-database --output-path {params.output_path} --use-parquet --parquet-dir {params.parquet_dir} --overwrite"
