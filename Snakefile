"""Snakemake pipeline for mitochondrial morphology screen analysis.

USAGE
=====
Three methods for generating virtual screen results:

Method 0: BASELINE (Recommended for production/publication)
  - Uses pre-computed validated results from S3
  - Fast, no recalculation
  - Run: snakemake all_baseline -c4 -p
  - Output: data/processed/screen_results_baseline.duckdb

Method 1: NOTEBOOK (Original exploratory implementation)
  - Recalculates from raw per-site profiles
  - Converted Jupyter notebook (retained for reference)
  - Run: snakemake all_notebook -c4 -p
  - Output: data/processed/screen_results_notebook.duckdb

Method 2: MODULE (Recommended for development)
  - Recalculates from raw per-site profiles
  - Clean refactored implementation
  - Includes diagnostic comparisons with baseline
  - Run: snakemake all_module -c4 -p
  - Output: data/processed/screen_results_module.duckdb

COMMON COMMANDS
===============
  snakemake --list-target-rules          # List all available targets
  snakemake all_baseline -c4 -p          # Method 0: Baseline pipeline
  snakemake all_module -c4 -p            # Method 2: Module pipeline
  snakemake all_module_diagnose -c4 -p   # Run diagnostics for all datasets
  snakemake -n -p <target>               # Dry run (preview)
  snakemake --summary                    # Show pipeline status

REPRODUCIBILITY
===============
Methods 1 & 2 regenerate results from raw data but show incomplete agreement with
baseline (cause unknown). Method 1 is the original implementation; Method 2 is a
refactored version with closer baseline agreement. Use Method 0 (baseline) for
validated publication results. See docs/PROGRESS.md for details.

CONFIGURATION
=============
Edit DATASETS variable below to process all 6 datasets or subset for testing.

"""

# ============================================================================
# Configuration
# ============================================================================

# Dataset configuration
# All 6 datasets available in the study
ALL_DATASETS = ["CDRP", "jump_compound", "jump_crispr", "jump_orf", "lincs", "taorf"]

# Filtered subset for Methods 1 & 2 (regeneration from raw data - slower)
# Method 0 (baseline) uses ALL_DATASETS since it only downloads pre-computed CSVs (fast)
DATASETS = ["lincs", "taorf"]

# S3 configuration
S3_BASE = "s3://imaging-platform/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace"

# Base paths (aligned with haghighi_mito/config.py)
DATA_DIR = "data"
EXTERNAL_DATA_DIR = f"{DATA_DIR}/external"
INTERIM_DATA_DIR = f"{DATA_DIR}/interim"
PROCESSED_DATA_DIR = f"{DATA_DIR}/processed"

# Mito project paths (from S3 download)
MITO_PROJECT_DIR = f"{EXTERNAL_DATA_DIR}/mito_project"
MITO_WORKSPACE_DIR = f"{MITO_PROJECT_DIR}/workspace"

# Method-specific output directories
BASELINE_DIR = f"{MITO_WORKSPACE_DIR}/results/virtual_screen_baseline"
NOTEBOOK_DIR = f"{PROCESSED_DATA_DIR}/virtual_screen_notebook"
MODULE_DIR = f"{PROCESSED_DATA_DIR}/virtual_screen_module"

# Intermediate Parquet directories
INTERIM_BASELINE = f"{INTERIM_DATA_DIR}/parquet_baseline"
INTERIM_NOTEBOOK = f"{INTERIM_DATA_DIR}/parquet_notebook"
INTERIM_MODULE = f"{INTERIM_DATA_DIR}/parquet_module"

# Processed tables directories
TABLES_BASELINE = f"{PROCESSED_DATA_DIR}/tables/generated_from_baseline"
TABLES_NOTEBOOK = f"{PROCESSED_DATA_DIR}/tables/generated_from_notebook"
TABLES_MODULE = f"{PROCESSED_DATA_DIR}/tables/generated_from_module"


# ============================================================================
# METHOD 0: BASELINE PIPELINE - Download + Process Pre-Computed Results
# ============================================================================
# This pipeline uses validated CSVs from July 2024 (uploaded to S3).
# It only performs filtering and formatting - NO slope calculation or stats.
# Output: data/processed/screen_results_baseline.duckdb (178,826 rows for all 6 datasets)
#
# Usage: snakemake all_baseline -c4 -p

## Download Rules ##

rule download_baseline_csv:
    """Download a single pre-computed virtual screen CSV from S3 (Method 0)."""
    output:
        csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    params:
        s3_path=lambda wildcards: f"{S3_BASE}/results/virtual_screen/{wildcards.dataset}_results_pattern_aug_070624.csv"
    shell:
        """
        mkdir -p {BASELINE_DIR}
        aws s3 cp {params.s3_path} {output.csv}
        """

## Processing Rules ##

rule process_baseline_csv:
    """Process baseline CSV to Excel + Parquet (filtering/formatting only)."""
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

rule create_baseline_database:
    """Combine all baseline Parquet files into unified DuckDB database."""
    input:
        expand(f"{INTERIM_BASELINE}/{{dataset}}_unfiltered.parquet",
               dataset=ALL_DATASETS)
    output:
        "data/processed/screen_results_baseline.duckdb"
    params:
        output_path="data/processed/screen_results_baseline.duckdb",
        datasets=",".join(ALL_DATASETS)
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_BASELINE} \
            --datasets {params.datasets} \
            --overwrite
        """

## Target Rules ##

rule all_baseline:
    """Target: Complete baseline pipeline (CSV → Excel → DuckDB)."""
    input:
        "data/processed/screen_results_baseline.duckdb"


# ============================================================================
# REGENERATED PIPELINES - Shared Data Downloads
# ============================================================================
# Both Method 1 (notebook) and Method 2 (clean module) require the same raw data:
# - Per-site aggregated profiles (2.7 GB across all 6 datasets)
# - Metadata files (~1.7 GB)
# - Orthogonal feature lists (~10 KB)
#
# Data is auto-downloaded when needed by Snakemake rules

rule download_orth_features:
    """Download all orthogonal feature lists (7 files, ~10 KB)."""
    output:
        touch(f"{MITO_WORKSPACE_DIR}/results/target_pattern_orth_features_lists/.download_complete")
    params:
        s3_dir=f"{S3_BASE}/results/target_pattern_orth_features_lists/",
        local_dir=f"{MITO_WORKSPACE_DIR}/results/target_pattern_orth_features_lists/"
    shell:
        """
        mkdir -p {params.local_dir}
        s5cmd sync '{params.s3_dir}*' {params.local_dir}
        """

rule download_per_site_profiles_dataset:
    """Download per-site profiles for a specific dataset (~100-2000 MB per dataset)."""
    output:
        touch(f"{MITO_WORKSPACE_DIR}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete")
    params:
        s3_dir=lambda wildcards: f"{S3_BASE}/per_site_aggregated_profiles_newpattern_2/{wildcards.dataset}/",
        local_dir=lambda wildcards: f"{MITO_WORKSPACE_DIR}/per_site_aggregated_profiles_newpattern_2/{wildcards.dataset}/"
    shell:
        """
        mkdir -p {params.local_dir}
        s5cmd sync '{params.s3_dir}*' {params.local_dir}
        """

rule download_metadata_file:
    """Download individual metadata file from S3."""
    output:
        f"{MITO_WORKSPACE_DIR}/metadata/{{metadata_path}}"
    params:
        s3_path=lambda wildcards: f"{S3_BASE}/metadata/{wildcards.metadata_path}"
    shell:
        """
        mkdir -p $(dirname {output})
        s5cmd cp {params.s3_path} {output}
        """

# ============================================================================
# METHOD 1: REGENERATED - Notebook (Original Implementation)
# ============================================================================
# Uses notebooks/2.0-mh-virtual-screen.py (converted Jupyter notebook).
# This is the original exploratory code, retained for reference.
#
# STATUS: ✅ COMPLETE - Full pipeline works (CSV → Excel → DuckDB)
#          Method 2 provides a cleaner refactored implementation
#
# Usage: snakemake all_notebook -c4 -p
# Output: data/processed/screen_results_notebook.duckdb

## Notebook Execution Rules ##

rule run_virtual_screen_notebook:
    """Run virtual screen using notebook 2.0 (recalculates slopes/stats)."""
    input:
        notebook="notebooks/2.0-mh-virtual-screen.py",
        # Ensure all required data is downloaded first
        orth_features=f"{MITO_WORKSPACE_DIR}/results/target_pattern_orth_features_lists/.download_complete",
        per_site_profiles=f"{MITO_WORKSPACE_DIR}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete",
        metadata_cdrp=f"{MITO_WORKSPACE_DIR}/metadata/CDRP_meta.csv",
        metadata_orf_list=f"{MITO_WORKSPACE_DIR}/metadata/JUMP-ORF/ORF_list.tsv",
        metadata_compound=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/compound.csv.gz",
        metadata_crispr=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/crispr.csv.gz",
        metadata_orf=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/orf.csv.gz",
        metadata_plate=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/plate.csv.gz",
        metadata_well=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/well.csv.gz",
        metadata_lincs=f"{MITO_WORKSPACE_DIR}/metadata/LINCS_meta.csv",
        metadata_taorf=f"{MITO_WORKSPACE_DIR}/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz",
        metadata_lincs_drug=f"{MITO_WORKSPACE_DIR}/metadata/lincs/DrugRepurposing_Metadata.csv"
    output:
        csv=f"{NOTEBOOK_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    shell:
        """
        mkdir -p {NOTEBOOK_DIR}
        pixi run python {input.notebook} --dataset {wildcards.dataset}
        """


## Processing Rules ##

rule process_notebook_csv:
    """Process notebook-generated CSV to Excel + Parquet (Method 1)."""
    input:
        csv=f"{NOTEBOOK_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        excel=f"{TABLES_NOTEBOOK}/{{dataset}}_screen_results.xlsx",
        parquet=f"{INTERIM_NOTEBOOK}/{{dataset}}_unfiltered.parquet"
    shell:
        """
        pixi run haghighi-mito process-csv-single \
            --dataset {wildcards.dataset} \
            --csv-path {input.csv} \
            --output-dir {TABLES_NOTEBOOK} \
            --parquet-output-dir {INTERIM_NOTEBOOK}
        """

rule create_notebook_database:
    """Combine notebook-generated Parquet files into unified DuckDB database (Method 1)."""
    input:
        expand(f"{INTERIM_NOTEBOOK}/{{dataset}}_unfiltered.parquet",
               dataset=DATASETS)
    output:
        "data/processed/screen_results_notebook.duckdb"
    params:
        output_path="data/processed/screen_results_notebook.duckdb",
        datasets=",".join(DATASETS)
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_NOTEBOOK} \
            --datasets {params.datasets} \
            --overwrite
        """

## Target Rules ##

rule all_notebook_csvs:
    """Target: Generate results CSVs for all datasets (notebook analysis only, no Excel/DuckDB)."""
    input:
        expand(f"{NOTEBOOK_DIR}/{{dataset}}_results_pattern_aug_070624.csv",
               dataset=DATASETS)

rule all_notebook:
    """Target: Complete Method 1 pipeline (notebook → CSV → Excel → DuckDB)."""
    input:
        "data/processed/screen_results_notebook.duckdb"


# ============================================================================
# METHOD 2: REGENERATED - Clean Module (Refactored Implementation)
# ============================================================================
# Uses haghighi_mito/virtual_screen.py (clean, documented implementation).
# This is a refactored version of the notebook logic.
#
# STATUS: ✅ COMPLETE - Full pipeline works (CSV → Excel → DuckDB)
#          Clean refactored implementation, recommended for development
#
# Usage: snakemake all_module -c4 -p
# Output: data/processed/screen_results_module.duckdb

## Analysis Rules ##

rule run_virtual_screen_module:
    """Run virtual screen using clean module (recalculates slopes/stats from per-site profiles)."""
    input:
        # Ensure all required data is downloaded first
        orth_features=f"{MITO_WORKSPACE_DIR}/results/target_pattern_orth_features_lists/.download_complete",
        per_site_profiles=f"{MITO_WORKSPACE_DIR}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete",
        metadata_cdrp=f"{MITO_WORKSPACE_DIR}/metadata/CDRP_meta.csv",
        metadata_orf_list=f"{MITO_WORKSPACE_DIR}/metadata/JUMP-ORF/ORF_list.tsv",
        metadata_compound=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/compound.csv.gz",
        metadata_crispr=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/crispr.csv.gz",
        metadata_orf=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/orf.csv.gz",
        metadata_plate=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/plate.csv.gz",
        metadata_well=f"{MITO_WORKSPACE_DIR}/metadata/JUMP/well.csv.gz",
        metadata_lincs=f"{MITO_WORKSPACE_DIR}/metadata/LINCS_meta.csv",
        metadata_taorf=f"{MITO_WORKSPACE_DIR}/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz",
        metadata_lincs_drug=f"{MITO_WORKSPACE_DIR}/metadata/lincs/DrugRepurposing_Metadata.csv"
    output:
        results_csv="data/processed/virtual_screen_module/{dataset}_results_pattern_aug_070624.csv"
    shell:
        """
        pixi run haghighi-mito virtual-screen --dataset {wildcards.dataset}
        """

rule diagnose_module:
    """Compare module CSV with baseline and generate diagnostic plots (FAST - ~1 sec)."""
    input:
        results_csv="data/processed/virtual_screen_module/{dataset}_results_pattern_aug_070624.csv",
        baseline_csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        comparison_csv="data/processed/virtual_screen_module/{dataset}_baseline_comparison.csv",
        plot="data/processed/virtual_screen_module/{dataset}_comparison_metrics.png"
    shell:
        """
        pixi run haghighi-mito compare-baseline --dataset {wildcards.dataset}
        """

## Processing Rules ##

rule process_module_csv:
    """Process module-generated CSV to Excel + Parquet (Method 2)."""
    input:
        csv=f"{MODULE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        excel=f"{TABLES_MODULE}/{{dataset}}_screen_results.xlsx",
        parquet=f"{INTERIM_MODULE}/{{dataset}}_unfiltered.parquet"
    shell:
        """
        pixi run haghighi-mito process-csv-single \
            --dataset {wildcards.dataset} \
            --csv-path {input.csv} \
            --output-dir {TABLES_MODULE} \
            --parquet-output-dir {INTERIM_MODULE}
        """

rule create_module_database:
    """Combine module-generated Parquet files into unified DuckDB database (Method 2)."""
    input:
        expand(f"{INTERIM_MODULE}/{{dataset}}_unfiltered.parquet",
               dataset=DATASETS)
    output:
        "data/processed/screen_results_module.duckdb"
    params:
        output_path="data/processed/screen_results_module.duckdb",
        datasets=",".join(DATASETS)
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_MODULE} \
            --datasets {params.datasets} \
            --overwrite
        """

## Target Rules ##

rule all_module_csvs:
    """Target: Generate results CSVs for all datasets (virtual screen analysis only, no Excel/DuckDB)."""
    input:
        expand("data/processed/virtual_screen_module/{dataset}_results_pattern_aug_070624.csv",
               dataset=DATASETS)

rule all_module_diagnose:
    """Target: Run diagnostics (comparison CSV + plots) for all datasets."""
    input:
        expand("data/processed/virtual_screen_module/{dataset}_comparison_metrics.png",
               dataset=DATASETS)

rule all_module:
    """Target: Complete Method 2 pipeline (CSV → Excel → DuckDB)."""
    input:
        "data/processed/screen_results_module.duckdb"


# ============================================================================
# UTILITY RULES
# ============================================================================

rule generate_dag_visualizations:
    """Generate DAG PNGs for all pipeline methods in docs/pipeline/."""
    shell:
        "./scripts/generate-dag.sh"


rule validate_database_pair:
    """Compare two DuckDB databases (e.g., baseline vs module).

    Usage: snakemake validate_database_pair --config db1=baseline db2=module
    """
    params:
        db1=config.get("db1", "baseline"),
        db2=config.get("db2", "module")
    shell:
        """
        pixi run haghighi-mito validate-databases \
            --baseline data/processed/screen_results_{params.db1}.duckdb \
            --new data/processed/screen_results_{params.db2}.duckdb
        """


rule run_reproduction_script:
    """Run slope discrepancy reproduction script for a dataset.

    Generates detailed comparison plots and CSV for debugging/author consultation.
    """
    input:
        baseline_csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv",
        per_site_profiles=f"{MITO_WORKSPACE_DIR}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete"
    output:
        plot="data/processed/virtual_screen_module/{dataset}_slope_discrepancy.png",
        csv="data/processed/virtual_screen_module/{dataset}_slope_comparison.csv"
    shell:
        """
        pixi run python scripts/reproduce_slope_discrepancy.py {wildcards.dataset}
        """


rule all_reproduce:
    """Target: Run reproduction script for all datasets."""
    input:
        expand("data/processed/virtual_screen_module/{dataset}_slope_discrepancy.png",
               dataset=ALL_DATASETS)


# ============================================================================
# Configuration Display
# ============================================================================

# Display configuration at pipeline start
onstart:
    print("=" * 70)
    print("MITOCHONDRIAL MORPHOLOGY SCREEN PIPELINE")
    print("=" * 70)
    print(f"Baseline datasets (Method 0): {', '.join(ALL_DATASETS)}")
    print(f"Module/Notebook datasets (Methods 1 & 2): {', '.join(DATASETS)}")
    print(f"Data directory: {DATA_DIR}")
    print(f"Mito workspace: {MITO_WORKSPACE_DIR}")
    print(f"Baseline output: {PROCESSED_DATA_DIR}/screen_results_baseline.duckdb")
    print(f"Notebook output: {PROCESSED_DATA_DIR}/screen_results_notebook.duckdb")
    print(f"Module output:   {PROCESSED_DATA_DIR}/screen_results_module.duckdb")
    print("=" * 70)
