"""Snakemake pipeline for mitochondrial morphology screen analysis."""

# ============================================================================
# Configuration
# ============================================================================

# Dataset configuration
DATASETS = ["CDRP", "jump_compound", "jump_crispr", "jump_orf", "lincs", "taorf"]

# S3 configuration
S3_BASE = "s3://imaging-platform/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace"

# Local directory structure
EXTERNAL_BASE = "data/external/mito_project/workspace"
BASELINE_DIR = f"{EXTERNAL_BASE}/results/virtual_screen_baseline"
REGEN_DIR = f"{EXTERNAL_BASE}/results/virtual_screen_regenerated"

# Processing output directories
INTERIM_BASELINE = "data/interim/parquet_baseline"
INTERIM_REGEN = "data/interim/parquet_regenerated"
TABLES_BASELINE = "data/processed/tables/generated_from_s3_baseline"
TABLES_REGEN = "data/processed/tables/generated_from_local"


# ============================================================================
# Targets
# ============================================================================

# Default rule - runs baseline pipeline
rule all:
    input:
        "data/processed/screen_results_baseline.duckdb"


# Target: Download all baseline CSVs from S3
rule download_all_baseline:
    input:
        expand(f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv",
               dataset=DATASETS)

# Target: Download all data needed for virtual screening analysis (notebook 2.0)
rule download_screening_data:
    input:
        f"{EXTERNAL_BASE}/results/target_pattern_orth_features_lists/.download_complete",
        expand(f"{EXTERNAL_BASE}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete",
               dataset=DATASETS),
        f"{EXTERNAL_BASE}/metadata/CDRP_meta.csv",
        f"{EXTERNAL_BASE}/metadata/JUMP-ORF/ORF_list.tsv",
        f"{EXTERNAL_BASE}/metadata/JUMP/compound.csv.gz",
        f"{EXTERNAL_BASE}/metadata/JUMP/crispr.csv.gz",
        f"{EXTERNAL_BASE}/metadata/JUMP/orf.csv.gz",
        f"{EXTERNAL_BASE}/metadata/JUMP/plate.csv.gz",
        f"{EXTERNAL_BASE}/metadata/JUMP/well.csv.gz",
        f"{EXTERNAL_BASE}/metadata/LINCS_meta.csv",
        f"{EXTERNAL_BASE}/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz",
        f"{EXTERNAL_BASE}/metadata/lincs/DrugRepurposing_Metadata.csv"


# ============================================================================
# Download Rules - Baseline Data from S3
# ============================================================================

# Download a single baseline CSV from S3
rule download_baseline_csv:
    output:
        csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    params:
        s3_path=lambda wildcards: f"{S3_BASE}/results/virtual_screen/{wildcards.dataset}_results_pattern_aug_070624.csv"
    shell:
        """
        mkdir -p {BASELINE_DIR}
        aws s3 cp {params.s3_path} {output.csv}
        """


# ============================================================================
# Download Rules - Virtual Screening Analysis Data (Notebook 2.0)
# ============================================================================

# Download all orthogonal feature lists (7 files, ~9.6 KB)
rule download_orth_features:
    output:
        touch(f"{EXTERNAL_BASE}/results/target_pattern_orth_features_lists/.download_complete")
    params:
        s3_dir=f"{S3_BASE}/results/target_pattern_orth_features_lists/",
        local_dir=f"{EXTERNAL_BASE}/results/target_pattern_orth_features_lists/"
    shell:
        """
        mkdir -p {params.local_dir}
        s5cmd sync '{params.s3_dir}*' {params.local_dir}
        """

# Download per-site profiles for a specific dataset
rule download_per_site_profiles_dataset:
    output:
        touch(f"{EXTERNAL_BASE}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete")
    params:
        s3_dir=lambda wildcards: f"{S3_BASE}/per_site_aggregated_profiles_newpattern_2/{wildcards.dataset}/",
        local_dir=lambda wildcards: f"{EXTERNAL_BASE}/per_site_aggregated_profiles_newpattern_2/{wildcards.dataset}/"
    shell:
        """
        mkdir -p {params.local_dir}
        s5cmd sync '{params.s3_dir}*' {params.local_dir}
        """

# Download individual metadata files
rule download_metadata_file:
    output:
        f"{EXTERNAL_BASE}/metadata/{{metadata_path}}"
    params:
        s3_path=lambda wildcards: f"{S3_BASE}/metadata/{wildcards.metadata_path}"
    shell:
        """
        mkdir -p $(dirname {output})
        s5cmd cp {params.s3_path} {output}
        """


# ============================================================================
# Processing Rules - Baseline Pipeline
# ============================================================================

# Process a single baseline CSV file to Excel and Parquet
rule process_baseline_csv:
    input:
        csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        excel=f"{TABLES_BASELINE}/{{dataset}}_screen_results.xlsx",
        parquet=f"{INTERIM_BASELINE}/{{dataset}}_unfiltered.parquet"
    shell:
        """
        pixi run haghighi-mito process-csv-single \
            --dataset {wildcards.dataset} \
            --csv-path {input.csv} \
            --output-dir {TABLES_BASELINE} \
            --parquet-output-dir {INTERIM_BASELINE}
        """


# Create unified DuckDB database from baseline Parquet files
rule create_baseline_database:
    input:
        expand(f"{INTERIM_BASELINE}/{{dataset}}_unfiltered.parquet",
               dataset=DATASETS)
    output:
        "data/processed/screen_results_baseline.duckdb"
    params:
        output_path="data/processed/screen_results_baseline.duckdb"
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_BASELINE} \
            --overwrite
        """


# ============================================================================
# Processing Rules - Regenerated Pipeline (for future use)
# ============================================================================

# Process a single regenerated CSV file to Excel and Parquet
rule process_regenerated_csv:
    input:
        csv=f"{REGEN_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        excel=f"{TABLES_REGEN}/{{dataset}}_screen_results.xlsx",
        parquet=f"{INTERIM_REGEN}/{{dataset}}_unfiltered.parquet"
    shell:
        """
        pixi run haghighi-mito process-csv-single \
            --dataset {wildcards.dataset} \
            --csv-path {input.csv} \
            --output-dir {TABLES_REGEN} \
            --parquet-output-dir {INTERIM_REGEN}
        """


# Create unified DuckDB database from regenerated Parquet files
rule create_regenerated_database:
    input:
        expand(f"{INTERIM_REGEN}/{{dataset}}_unfiltered.parquet",
               dataset=DATASETS)
    output:
        "data/processed/screen_results_regenerated.duckdb"
    params:
        output_path="data/processed/screen_results_regenerated.duckdb"
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_REGEN} \
            --overwrite
        """


# ============================================================================
# Additional Targets
# ============================================================================

# Target: Run regenerated pipeline
rule all_regenerated:
    input:
        "data/processed/screen_results_regenerated.duckdb"


# ============================================================================
# Configuration Display
# ============================================================================

# Display configuration at pipeline start
onstart:
    print("=" * 70)
    print("MITOCHONDRIAL MORPHOLOGY SCREEN PIPELINE")
    print("=" * 70)
    print(f"Datasets: {', '.join(DATASETS)}")
    print(f"Baseline input:  {BASELINE_DIR}")
    print(f"Baseline output: {TABLES_BASELINE}")
    print(f"Baseline Parquet: {INTERIM_BASELINE}")
    print(f"Baseline DB: data/processed/screen_results_baseline.duckdb")
    print(f"Regenerated input:  {REGEN_DIR}")
    print(f"Regenerated output: {TABLES_REGEN}")
    print(f"Regenerated Parquet: {INTERIM_REGEN}")
    print(f"Regenerated DB: data/processed/screen_results_regenerated.duckdb")
    print("=" * 70)
