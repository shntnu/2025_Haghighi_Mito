# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis pipelines and results for a mitochondrial morphology study using Cell Painting assays. The project performs virtual screens across multiple perturbation datasets (compounds, ORFs, CRISPR) to identify genetic and chemical modulators of mitochondrial radial distribution patterns in fibroblasts.

## Repository Structure

This repository follows the [Carpenter-Singh lab workflow conventions](protocols/workflows.md) based on Cookiecutter Data Science principles:

```
2025_Haghighi_Mito/
├── data/                  # Project data following lab conventions
│   ├── external/          # Downloaded data from S3 (gitignored)
│   ├── interim/           # Intermediate processing outputs (Parquet files)
│   └── processed/         # Final analysis outputs
│       ├── figures/       # Publication figures (PDFs)
│       └── tables/        # Screen results (XLSX) - versioned directories
├── haghighi_mito/         # Python package with processing code
├── notebooks/             # Analysis notebooks (numbered by phase)
├── pipelines/             # CellProfiler pipelines (.cppipe files)
├── protocols/             # Experimental protocols and workflow documentation
├── docs/                  # Documentation and analysis notes
│   ├── DATA_DOWNLOAD.md   # S3 data download requirements
│   ├── DATA_FLOW.md       # Complete pipeline data flow documentation
│   └── PROGRESS.md        # High-level progress log
├── Snakefile              # Pipeline automation (download + processing)
├── Justfile               # Convenience commands (run, clean, etc.)
├── manuscript.md          # Manuscript content
├── CLAUDE.md             # Repository documentation
└── README.md             # Project overview
```

### Workflow Principles

Following the conventions documented in `protocols/workflows.md`:

1. **Data flows in one direction**: `raw/` + `external/` → `interim/` → `processed/`
2. **Raw data is immutable**: Never edit raw data directly
3. **Clear separation of concerns**: Pipeline-managed data vs. personal analysis
4. **Everything is reproducible**: Code + raw data = any output

### Analysis Notebooks (`notebooks/`)

Notebooks follow Carpenter-Singh lab convention: `<phase>.<sequence>-<initials>-<description>.py`

- Phase 1 = Data cleaning/feature engineering
- Phase 2 = Analysis/screening

- **`1.0-mh-feat-importance.py`**: Feature importance analysis for mitochondrial morphology measurements.

- **`2.0-mh-virtual-screen.py`**: Primary pipeline for screening phenotype strength of mitochondrial radial distribution features across datasets (LINCS, CDRP, JUMP-ORF, JUMP-CRISPR, JUMP-Compound, TA-ORF). Performs statistical testing, filtering, and generates ranked perturbation lists.

- **`2.1-mh-set-enrichment-analysis.py`**: Gene set enrichment analysis (GSEA) using blitzgsea package. Analyzes enrichment in GO terms, pathways (KEGG, WikiPathways), disease associations (OMIM, HPO), and mechanism-of-action (MOA) categories.

- **`2.2-mh-check-vs-lists.py`**: Validation and comparison of virtual screen hit lists.

### CellProfiler Pipelines (`pipelines/`)

All pipelines are `.cppipe` files for CellProfiler image analysis:

- **`dapi_mito_actin_fibroblasts_40x/`**: Pipelines for reprocessing the repurposing dataset
  - `illum.cppipe` / `illum_without_batchfile.cppipe`: Illumination correction
  - `analysis.cppipe` / `analysis_without_batchfile.cppipe`: Feature extraction

- **Main analysis pipelines**:
  - `analysis_with_images*.cppipe`: Analysis with image export (various versions for CP3/batch processing)
  - `analysis_bethac07_radial_entropy*.cppipe`: Radial entropy-focused analysis
  - `radialdistribution*.cppipe`: Radial distribution measurements
  - `qc*.cppipe`: Quality control pipelines

- **Cell-type specific**:
  - `neurons_segmentation_CPA.cppipe`: Neuron segmentation
  - `oligodendrocytes.cppipe`: Oligodendrocyte analysis

### Data (`data/`)

Organized following Carpenter-Singh lab data flow conventions:

- **`data/external/mito_project/workspace/`**: Downloaded data from S3 (gitignored)
  - `results/virtual_screen_baseline/`: Pre-computed virtual screen CSVs (July 2024, validated)
  - `results/virtual_screen_regenerated/`: Locally-generated CSVs from notebook 2.0 (experimental)
  - `per_site_aggregated_profiles_newpattern_2/`: Per-site profile data by dataset
  - `metadata/`: Dataset metadata files
  - Reference datasets: `KEGG_2021_Human_table.txt`, `WikiPathways_2024_Human_table.txt`

- **`data/interim/`**: Intermediate pipeline outputs (gitignored)
  - `parquet_baseline/`: Parquet files from baseline CSVs
  - `parquet_regenerated/`: Parquet files from regenerated CSVs

- **`data/processed/`**: Final analysis outputs
  - **`data/processed/tables/`**: Virtual screen results - versioned directories
    - `curated_2024-08-11/`: Original curated Excel files (git-tracked)
    - `curated_2025-10-25/`: Current curated Excel files (git-tracked)
    - `generated_from_s3_baseline/`: Pipeline outputs from baseline CSVs (gitignored, reproducible)
    - `generated_from_local/`: Pipeline outputs from regenerated CSVs (gitignored, experimental)
  - **`data/processed/figures/`**: Publication figures (PDFs)
  - **`screen_results_baseline.duckdb`**: Query interface for baseline results (178,826 rows)

### Protocols (`protocols/`)

- **`McleanCollectionFibroblastGrowthProtocol.md`**: Detailed protocol for fibroblast culture, passaging, staining, and imaging

## Key Analysis Features

### Virtual Screen Pipeline

The virtual screen identifies perturbations affecting mitochondrial radial distribution:

1. **Target Feature**: Primary metric is `slope` (radial distribution slope) or `d_slope` (effect size)
2. **Statistical Testing**: Per-site aggregation followed by per-plate statistical tests vs. controls
3. **Filtering Strategy**:
   - Cell count filter: Remove bottom 10% by cell count
   - Target feature significance: Benjamini-Hochberg corrected p-values
   - Orthogonal feature filter: Ensure perturbation doesn't affect unrelated features
4. **Datasets Analyzed**: LINCS (compounds), CDRP (compounds), JUMP-ORF/CRISPR (genetics), TA-ORF (genetics)

### Dataset Parameters

Each dataset has specific metadata columns and identifiers defined in `datasets_info_dict`:

- **LINCS**: Key column `Metadata_pert_name`, reference set `Metadata_pert_id_dose`
- **CDRP**: Key column `Metadata_pert_id`, reference set `Metadata_Sample_Dose`
- **JUMP-ORF**: Key column `Metadata_Symbol`, reference set `Metadata_JCP2022`
- **JUMP-CRISPR**: Key column `Metadata_Symbol`, reference set `Metadata_JCP2022`
- **JUMP-Compound**: Key column `Metadata_InChIKey`, reference set `Metadata_JCP2022`
- **TA-ORF**: Key column `Metadata_gene_name`, reference set `Metadata_broad_sample`

### Enrichment Analysis

Set enrichment uses **blitzgsea** package for:

- Gene ontology (GO Biological Process, Molecular Function, Cellular Component)
- Pathways (KEGG, WikiPathways, PFOCR)
- Disease associations (OMIM, HPO, GWAS Catalog)
- MOA categories for compound screens
- Custom gene lists (mitochondrial dynamics, fusion/fission/transport/mitophagy)

## Working with This Repository

### Quick Start - Running the Pipeline

The repository uses **Snakemake** for pipeline automation and **Justfile** for convenience commands.

**View available commands:**
```bash
just --list
```

**Process baseline data (recommended - validated results):**
```bash
just download-baseline  # Download pre-computed virtual screen CSVs from S3
just run-baseline       # Process CSVs → Excel + Parquet → DuckDB
```

**Download raw data for notebook 2.0:**
```bash
just download-raw       # Download metadata, profiles, orth features (~2.8 GB)
```

**Visualize pipeline DAG:**
```bash
just viz                # Generate pipeline diagrams in docs/pipeline/
```

### Python Package (`haghighi_mito/`)

Processing code is organized as a package with CLI interface:

```bash
pixi run haghighi-mito --help              # View all commands
pixi run haghighi-mito process-csv-single  # Process single dataset CSV
pixi run haghighi-mito create-database     # Create DuckDB from Parquet files
pixi run haghighi-mito validate-databases  # Compare two DuckDB databases
```

Package modules:
- `config.py`: Dataset configuration (DATASET_INFO dict)
- `data.py`: Processing functions (CSV → Excel/Parquet, DuckDB creation)
- `vectorized_slope.py`: Optimized slope calculation (~200x speedup)
- `vectorized_stats.py`: Batch statistical testing
- `cli.py`: Typer CLI interface

### Analysis Notebooks

Notebooks are converted to `.py` scripts for command-line execution:

- **`1.0-mh-feat-importance.py`**: Feature importance analysis
- **`2.0-mh-virtual-screen.py`**: Virtual screening pipeline (accepts `--dataset` argument)
- **`2.1-mh-set-enrichment-analysis.py`**: Gene set enrichment analysis (GSEA)
- **`2.2-mh-check-vs-lists.py`**: Validation and comparison

### Key Dependencies

- **Pipeline**: snakemake, pixi, s5cmd (fast S3 downloads)
- **Data**: pandas, numpy, scipy, pyarrow (Parquet), duckdb
- **Visualization**: matplotlib, seaborn
- **Analysis**: sklearn, blitzgsea
- **Custom**: singlecell package (morphological analysis)

## Critical Analysis Parameters

- **BH-corrected critical values** (FDR=0.05): Dataset-specific, stored in `orth_bh_corrected_critical_dict` and `target_bh_corrected_critical_dict`
- **Minimum gene set size**: 4-10 depending on analysis
- **Cell count quantile filter**: 0.1 (bottom 10%)
- **Target feature**: `slope` or `d_slope` (effect size of radial distribution slope)

## Known Issues

### Baseline Non-Reproducibility

**Problem:** Locally-regenerated results differ from S3 baseline (July 2024) by ~77% (within 10% tolerance).

**Root Cause:** Baseline was generated with code that **no longer exists** in this repository:
- Repository created September 2025, baseline uploaded July 2024
- Current code (both `if 0` and `if 1` branches) produces different results
- 100% of large differences have different `last_peak_ind` values (peak detection differs)
- Input data confirmed identical (Count_Cells_avg matches 100%)
- Vectorization confirmed correct (produces identical results to row-by-row implementation)

**Impact:**
- Baseline pipeline (`just run-baseline`) uses validated S3 CSVs → **safe for publication**
- Regenerated pipeline produces experimental results → **use for validation only**
- Optimization work (vectorization) cannot be validated against baseline

**Decision Point:** Accept baseline as-is OR commit to regenerated version (requires domain expert validation). See docs/PROGRESS.md for detailed analysis.

**Recommendation:** Always use `just run-baseline` for validated results. Regenerated pipeline exists for future improvements.

### Data Provenance

Excel files in `curated_2024-08-11/` and `curated_2025-10-25/` are **manually curated via Google Sheets**, not direct pipeline outputs. Pipeline-generated Excel files are in `generated_from_*/` directories (gitignored, reproducible).

## Important Instructions for Claude Code

- Do what has been asked; nothing more, nothing less
- NEVER create files unless absolutely necessary for achieving your goal
- ALWAYS prefer editing an existing file to creating a new one
- NEVER proactively create documentation files (*.md) or README files unless explicitly requested
- NEVER use emojis in any output, documentation, code comments, or commit messages unless explicitly requested
- Always use `pixi run python` instead of just `python`
- DO NOT update docs/PROGRESS.md unless explicitly asked by the user (see maintenance guidelines in that file)