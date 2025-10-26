"""Snakemake pipeline for mitochondrial morphology screen analysis.

================================================================================
PIPELINE OVERVIEW: Three Methods for Virtual Screen Results
================================================================================

This pipeline supports THREE methods for generating virtual screen results.
Each method produces different outputs with different trade-offs.

METHOD 0: BASELINE (July 2024 - Validated, Fast)
────────────────────────────────────────────────
What:     Pre-computed CSVs uploaded to S3 in July 2024 (65 MB total)
Status:   ✅ COMPLETE - Validated in published manuscript
Use case: Production analysis, publication results
Code:     haghighi_mito/data.py (formatting/filtering only - NO recalculation)
Commands: just download-baseline && just run-baseline
Output:   data/processed/screen_results_baseline.duckdb (178,826 rows)
Speed:    ~5 minutes

Data flow:
  S3 CSVs (slopes already computed)
    → Download (aws s3 cp)
    → Filter + format to Excel/Parquet (data.py::process_single_virtual_screen_csv)
    → Create DuckDB (data.py::create_screen_database)

METHOD 1: REGENERATED - Notebook (Complete but Messy)
──────────────────────────────────────────────────────
What:     Exploratory Jupyter notebook converted to .py script (1,433 lines)
Status:   ✅ COMPLETE - Full pipeline works (CSV → Excel → DuckDB)
Use case: Only regenerated method with complete output chain
Code:     notebooks/2.0-mh-virtual-screen.py + haghighi_mito/data.py
Commands: just download-raw && just run-regenerated
Output:   data/processed/screen_results_regenerated.duckdb
Speed:    ~10 minutes per dataset
Note:     Messy code with dead branches (if 0:), but functional
          NEEDED until Method 2 gains Excel/DuckDB processing steps

Data flow:
  Raw per-site profiles (2.7 GB, 17K rows for taorf)
    → Download (s5cmd sync)
    → Calculate slopes, peaks, stats (notebook 2.0)
    → Save CSV (328 rows for taorf)
    → Filter + format to Excel/Parquet (data.py::process_single_virtual_screen_csv)
    → Create DuckDB (data.py::create_screen_database)

METHOD 2: REGENERATED - Clean Module (Active Development, Incomplete)
──────────────────────────────────────────────────────────────────────
What:     Professional refactor of notebook logic (448 clean lines)
Status:   ⚠️ INCOMPLETE - Stops at CSV generation, missing Excel/DuckDB steps
Gap:      Needs process_csv + create_database rules (exist for Method 1)
Use case: Baseline comparison, diagnostics, methodology validation
Code:     haghighi_mito/virtual_screen.py + vectorized helpers + diagnostics.py
Commands: just download-raw && just run-virtual-screen-for DATASET
Output:   CSVs in virtual_screen_simple/ + baseline comparison + plots
Speed:    ~10 minutes per dataset
TODO:     Add Excel/Parquet/DuckDB processing → then can deprecate Method 1

Data flow (current):
  Raw per-site profiles (2.7 GB, 17K rows for taorf)
    → Download (s5cmd sync)
    → Calculate slopes, peaks, stats (virtual_screen.py::run_virtual_screen)
    → Compare to baseline (diagnostics.py::compare_with_baseline)
    → Save CSV + comparison + plots
    → [STOPS HERE - missing Excel/DuckDB steps]

Data flow (future):
  ... (same as above) →
    → Filter + format to Excel/Parquet (data.py::process_single_virtual_screen_csv)
    → Create DuckDB (data.py::create_screen_database)

================================================================================
THE REPRODUCIBILITY ISSUE
================================================================================

All regenerated methods produce ~77% agreement with baseline due to lost
July 2024 code. Key findings:
- Baseline generated with code that predates this repository (Sept 2025)
- Input data confirmed identical (Count_Cells_avg matches 100%)
- 100% of large slope differences have different peak indices
- Root cause: Unknown peak detection algorithm used in July 2024
- Impact: Baseline trusted but opaque; regenerated methods transparent but inexact

Match rates (taorf, n=327):
- t_target_pattern: 37% within 10% (no peak detection, best metric)
- slope: 19% within 10% (peak detection differs)
- t_orth: 33% within 10% (validates processing)

See docs/PROGRESS.md for detailed investigation history (Oct 2025).

================================================================================
OUTPUT DIRECTORY STRUCTURE
================================================================================

data/
├── external/mito_project/workspace/results/
│   ├── virtual_screen_baseline/          # Method 0: S3 downloads (65 MB)
│   └── virtual_screen_regenerated/       # Method 1: Notebook CSVs
│
├── processed/
│   ├── screen_results_baseline.duckdb    # Method 0 final database
│   ├── screen_results_regenerated.duckdb # Method 1 final database
│   ├── screen_results_module.duckdb      # Method 2 final database (TODO)
│   │
│   ├── tables/
│   │   ├── generated_from_s3_baseline/   # Method 0 Excel files
│   │   ├── generated_from_local/         # Method 1 Excel files
│   │   └── generated_from_module/        # Method 2 Excel files (TODO)
│   │
│   ├── virtual_screen_simple/            # Method 2: CSVs + comparisons
│   │   ├── {dataset}_results_pattern_aug_070624.csv
│   │   └── {dataset}_baseline_comparison.csv
│   │
│   └── figures/t_target_pattern_analysis/  # Method 2: Diagnostic plots
│       └── {dataset}_baseline_vs_regenerated.png
│
└── interim/
    ├── parquet_baseline/                 # Method 0 intermediate
    ├── parquet_regenerated/              # Method 1 intermediate
    └── parquet_module/                   # Method 2 intermediate (TODO)

"""

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
NOTEBOOK_DIR = f"{EXTERNAL_BASE}/results/virtual_screen_regenerated"  # Method 1 output
MODULE_DIR = "data/processed/virtual_screen_simple"  # Method 2 output

# Processing output directories
INTERIM_BASELINE = "data/interim/parquet_baseline"
INTERIM_NOTEBOOK = "data/interim/parquet_regenerated"  # Method 1 intermediate
INTERIM_MODULE = "data/interim/parquet_module"  # Method 2 intermediate
TABLES_BASELINE = "data/processed/tables/generated_from_s3_baseline"
TABLES_NOTEBOOK = "data/processed/tables/generated_from_local"  # Method 1 tables
TABLES_MODULE = "data/processed/tables/generated_from_module"  # Method 2 tables


# ============================================================================
# Default Targets
# ============================================================================

# Default rule - runs baseline pipeline (Method 0)
rule all_baseline:
    input:
        "data/processed/screen_results_baseline.duckdb"

# Alias for default target
rule all:
    input:
        rules.all_baseline.input


# ============================================================================
# METHOD 0: BASELINE PIPELINE - Download + Process Pre-Computed Results
# ============================================================================
# This pipeline uses validated CSVs from July 2024 (uploaded to S3).
# It only performs filtering and formatting - NO slope calculation or stats.
# Output: data/processed/screen_results_baseline.duckdb (178,826 rows)
#
# Commands: just download-baseline && just run-baseline

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

rule download_all_baseline:
    """Target: Download all 6 baseline CSVs from S3 (65 MB total)."""
    input:
        expand(f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv",
               dataset=DATASETS)

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
# REGENERATED PIPELINES - Shared Data Downloads
# ============================================================================
# Both Method 1 (notebook) and Method 2 (clean module) require the same raw data:
# - Per-site aggregated profiles (2.7 GB across 6 datasets)
# - Metadata files (~1.7 GB)
# - Orthogonal feature lists (~10 KB)
#
# Commands: just download-raw

rule download_orth_features:
    """Download all orthogonal feature lists (7 files, ~10 KB)."""
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

rule download_per_site_profiles_dataset:
    """Download per-site profiles for a specific dataset (~100-2000 MB per dataset)."""
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

rule download_metadata_file:
    """Download individual metadata file from S3."""
    output:
        f"{EXTERNAL_BASE}/metadata/{{metadata_path}}"
    params:
        s3_path=lambda wildcards: f"{S3_BASE}/metadata/{wildcards.metadata_path}"
    shell:
        """
        mkdir -p $(dirname {output})
        s5cmd cp {params.s3_path} {output}
        """

rule download_screening_data:
    """Target: Download all raw data for virtual screening (2.7 GB total)."""
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
# METHOD 1: REGENERATED - Notebook (Complete but Messy)
# ============================================================================
# Uses notebooks/2.0-mh-virtual-screen.py (1,433 lines of converted Jupyter).
# This is the original exploratory code with lots of dead branches (if 0:).
#
# STATUS: ✅ COMPLETE - Full pipeline works (CSV → Excel → DuckDB)
#          NEEDED until Method 2 gains Excel/DuckDB processing steps
#
# Commands: just download-raw && just run-regenerated
# Output: data/processed/screen_results_regenerated.duckdb

## Notebook Execution Rules ##

rule run_virtual_screen_notebook:
    """Run virtual screen using notebook 2.0 (1433 lines, recalculates slopes/stats)."""
    input:
        notebook="notebooks/2.0-mh-virtual-screen.py",
        # Ensure all required data is downloaded first
        orth_features=f"{EXTERNAL_BASE}/results/target_pattern_orth_features_lists/.download_complete",
        per_site_profiles=f"{EXTERNAL_BASE}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete",
        metadata_cdrp=f"{EXTERNAL_BASE}/metadata/CDRP_meta.csv",
        metadata_orf_list=f"{EXTERNAL_BASE}/metadata/JUMP-ORF/ORF_list.tsv",
        metadata_compound=f"{EXTERNAL_BASE}/metadata/JUMP/compound.csv.gz",
        metadata_crispr=f"{EXTERNAL_BASE}/metadata/JUMP/crispr.csv.gz",
        metadata_orf=f"{EXTERNAL_BASE}/metadata/JUMP/orf.csv.gz",
        metadata_plate=f"{EXTERNAL_BASE}/metadata/JUMP/plate.csv.gz",
        metadata_well=f"{EXTERNAL_BASE}/metadata/JUMP/well.csv.gz",
        metadata_lincs=f"{EXTERNAL_BASE}/metadata/LINCS_meta.csv",
        metadata_taorf=f"{EXTERNAL_BASE}/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz",
        metadata_lincs_drug=f"{EXTERNAL_BASE}/metadata/lincs/DrugRepurposing_Metadata.csv"
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
        "data/processed/screen_results_regenerated.duckdb"
    params:
        output_path="data/processed/screen_results_regenerated.duckdb"
    shell:
        """
        pixi run haghighi-mito create-database \
            --output-path {params.output_path} \
            --use-parquet \
            --parquet-dir {INTERIM_NOTEBOOK} \
            --overwrite
        """

## Target Rules ##

rule all_notebook:
    """Target: Complete Method 1 pipeline (notebook → CSV → Excel → DuckDB)."""
    input:
        "data/processed/screen_results_regenerated.duckdb"

rule run_all_virtual_screen_notebooks:
    """Target: Run notebook method for all datasets (CSV generation only)."""
    input:
        expand(f"{NOTEBOOK_DIR}/{{dataset}}_results_pattern_aug_070624.csv",
               dataset=DATASETS)


# ============================================================================
# METHOD 2: REGENERATED - Clean Module (Active Development, Incomplete)
# ============================================================================
# Uses haghighi_mito/virtual_screen.py (448 clean lines of documented code).
# This is a professional refactor of the notebook logic.
#
# STATUS: ⚠️ INCOMPLETE - Stops at CSV generation + baseline comparison
#          Missing Excel + Parquet + DuckDB processing steps
#
# TODO: Add process_virtual_screen_simple_csv and create_virtual_screen_simple_database
#       rules (similar to Method 1) to complete the pipeline.
#       Once added, Method 1 (notebook) can be deprecated.
#
# Commands: just download-raw && just run-virtual-screen-for DATASET
# Current output: CSVs in virtual_screen_simple/ + comparison + plots

## Analysis Rules ##

rule run_virtual_screen_analysis:
    """Run virtual screen using clean module (448 lines, recalculates slopes/stats)."""
    input:
        # Ensure all required data is downloaded first
        orth_features=f"{EXTERNAL_BASE}/results/target_pattern_orth_features_lists/.download_complete",
        per_site_profiles=f"{EXTERNAL_BASE}/per_site_aggregated_profiles_newpattern_2/{{dataset}}/.download_complete",
        metadata_cdrp=f"{EXTERNAL_BASE}/metadata/CDRP_meta.csv",
        metadata_orf_list=f"{EXTERNAL_BASE}/metadata/JUMP-ORF/ORF_list.tsv",
        metadata_compound=f"{EXTERNAL_BASE}/metadata/JUMP/compound.csv.gz",
        metadata_crispr=f"{EXTERNAL_BASE}/metadata/JUMP/crispr.csv.gz",
        metadata_orf=f"{EXTERNAL_BASE}/metadata/JUMP/orf.csv.gz",
        metadata_plate=f"{EXTERNAL_BASE}/metadata/JUMP/plate.csv.gz",
        metadata_well=f"{EXTERNAL_BASE}/metadata/JUMP/well.csv.gz",
        metadata_lincs=f"{EXTERNAL_BASE}/metadata/LINCS_meta.csv",
        metadata_taorf=f"{EXTERNAL_BASE}/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz",
        metadata_lincs_drug=f"{EXTERNAL_BASE}/metadata/lincs/DrugRepurposing_Metadata.csv",
        # Also need baseline CSV for comparison
        baseline_csv=f"{BASELINE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
    output:
        results_csv="data/processed/virtual_screen_simple/{dataset}_results_pattern_aug_070624.csv",
        comparison_csv="data/processed/virtual_screen_simple/{dataset}_baseline_comparison.csv"
    shell:
        """
        pixi run haghighi-mito virtual-screen --dataset {wildcards.dataset} --compare-baseline
        """


rule plot_baseline_comparison:
    """Generate 2x2 scatter plots comparing regenerated vs baseline metrics."""
    input:
        comparison_csv="data/processed/virtual_screen_simple/{dataset}_baseline_comparison.csv"
    output:
        plot="data/processed/figures/t_target_pattern_analysis/{dataset}_baseline_vs_regenerated.png"
    shell:
        """
        pixi run haghighi-mito plot-baseline-comparison --dataset {wildcards.dataset}
        """


## Target Rules ##

rule run_all_virtual_screen_analysis:
    """Target: Run clean module method for all datasets (CSV generation + comparison)."""
    input:
        expand("data/processed/virtual_screen_simple/{dataset}_results_pattern_aug_070624.csv",
               dataset=DATASETS)

rule plot_all_baseline_comparisons:
    """Target: Generate comparison plots for all datasets."""
    input:
        expand("data/processed/figures/t_target_pattern_analysis/{dataset}_baseline_vs_regenerated.png",
               dataset=DATASETS)

## TODO: Missing Processing Rules for Method 2 Completion ##
#
# Uncomment and implement these rules to complete Method 2 pipeline.
# These rules mirror the Method 1 (notebook) pattern.
#
# rule process_module_csv:
#     """Process module-generated CSV to Excel + Parquet (Method 2)."""
#     input:
#         csv=f"{MODULE_DIR}/{{dataset}}_results_pattern_aug_070624.csv"
#     output:
#         excel=f"{TABLES_MODULE}/{{dataset}}_screen_results.xlsx",
#         parquet=f"{INTERIM_MODULE}/{{dataset}}_unfiltered.parquet"
#     shell:
#         """
#         pixi run haghighi-mito process-csv-single \
#             --dataset {wildcards.dataset} \
#             --csv-path {input.csv} \
#             --output-dir {TABLES_MODULE} \
#             --parquet-output-dir {INTERIM_MODULE}
#         """
#
# rule create_module_database:
#     """Combine module-generated Parquet files into unified DuckDB database (Method 2)."""
#     input:
#         expand(f"{INTERIM_MODULE}/{{dataset}}_unfiltered.parquet",
#                dataset=DATASETS)
#     output:
#         "data/processed/screen_results_module.duckdb"
#     params:
#         output_path="data/processed/screen_results_module.duckdb"
#     shell:
#         """
#         pixi run haghighi-mito create-database \
#             --output-path {params.output_path} \
#             --use-parquet \
#             --parquet-dir {INTERIM_MODULE} \
#             --overwrite
#         """
#
# rule all_module:
#     """Target: Complete Method 2 pipeline (CSV → Excel → DuckDB).
#
#     This is the Method 2 equivalent of all_notebook (Method 1).
#     Once implemented, use: just run-module
#     """
#     input:
#         "data/processed/screen_results_module.duckdb"


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
    print(f"Notebook input:  {NOTEBOOK_DIR}")
    print(f"Notebook output: {TABLES_NOTEBOOK}")
    print(f"Notebook Parquet: {INTERIM_NOTEBOOK}")
    print(f"Notebook DB: data/processed/screen_results_regenerated.duckdb")
    print("=" * 70)
