# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis pipelines and results for a mitochondrial morphology study using Cell Painting assays. The project performs virtual screens across multiple perturbation datasets (compounds, ORFs, CRISPR) to identify genetic and chemical modulators of mitochondrial radial distribution patterns in fibroblasts.

Supporting repository for:
> Haghighi, M. et al. Identifying and targeting abnormal mitochondrial localization associated with psychoses. *bioRxiv* 2025.10.08.676630 (2025) doi:[10.1101/2025.10.08.676630](https://doi.org/10.1101/2025.10.08.676630).

## Repository Structure

This repository follows the [Carpenter-Singh lab workflow conventions](protocols/workflows.md):

```
2025_Haghighi_Mito/
├── haghighi_mito/         # Python package with processing code
│   ├── cli.py             # Typer CLI interface
│   ├── config.py          # Dataset configuration (DATASET_INFO dict)
│   ├── data.py            # Processing: CSV → Excel/Parquet, DuckDB creation
│   ├── virtual_screen.py  # Clean module: per-site profiles → virtual screen CSVs
│   ├── diagnostics.py     # Baseline comparison plots and metrics
│   ├── vectorized_slope.py # Optimized slope calculation (~200x speedup)
│   └── vectorized_stats.py # Batch statistical testing
├── notebooks/             # Analysis notebooks (numbered by phase)
│   ├── 1.0-mh-feat-importance.py           # Feature importance
│   ├── 2.0-mh-virtual-screen.py            # Virtual screening (1433 lines)
│   ├── 2.1-mh-set-enrichment-analysis.py   # GSEA
│   └── 2.2-mh-check-vs-lists.py            # Validation
├── data/
│   ├── external/          # Downloaded data from S3 (gitignored)
│   ├── interim/           # Intermediate Parquet files (gitignored)
│   └── processed/         # Final outputs (DuckDB, Excel, figures)
├── Snakefile              # Pipeline automation (download + processing)
├── Justfile               # Convenience commands (run, clean, viz)
├── pyproject.toml         # Pixi configuration and dependencies
└── protocols/             # Experimental protocols and workflow docs
```

## Common Commands

### Running the Pipeline

View available commands:
```bash
just --list
```

**Recommended: Process validated baseline data (July 2024):**
```bash
just download-baseline  # Download pre-computed CSVs (65 MB, ~1 min)
just run-baseline       # Process → Excel + DuckDB (~5 min)
```

**Alternative: Regenerate from raw data (Method 2 - clean module):**
```bash
just download-raw           # Download per-site profiles (2.7 GB, ~5 min)
just run-module             # Run full pipeline → CSV → Excel → DuckDB (~10 min/dataset)
just run-module-for taorf   # Run single dataset
```

**Compare regenerated results with baseline:**
```bash
just compare-baseline-for taorf    # Compare CSVs (~1 sec)
just plot-comparison-for taorf     # Generate diagnostic plots
```

### Python CLI Commands

The `haghighi-mito` CLI provides lower-level control:

```bash
pixi run haghighi-mito --help                    # View all commands

# Run virtual screen from per-site profiles
pixi run haghighi-mito virtual-screen --dataset taorf

# Compare with baseline (fast, no regeneration)
pixi run haghighi-mito compare-baseline --dataset taorf

# Generate diagnostic plots
pixi run haghighi-mito plot-baseline-comparison --dataset taorf

# Process CSV to Excel + Parquet
pixi run haghighi-mito process-csv-single \
    --dataset taorf \
    --csv-path data/processed/virtual_screen_module/taorf_results_pattern_aug_070624.csv \
    --output-dir data/processed/tables/generated_from_module \
    --parquet-output-dir data/interim/parquet_module

# Create DuckDB from Parquet files
pixi run haghighi-mito create-database \
    --output-path data/processed/screen_results_module.duckdb \
    --use-parquet \
    --parquet-dir data/interim/parquet_module \
    --datasets lincs,taorf \
    --overwrite
```

### Development

```bash
pixi run python           # Use this instead of bare 'python'
pixi run jupyter lab      # Launch Jupyter
pixi run pytest           # Run tests
pixi run ruff check .     # Lint
pixi run ruff format .    # Format
```

## Pipeline Architecture

The repository supports **three methods** for generating virtual screen results:

### Method 0: Baseline (Validated, July 2024)
- **Status**: Production-ready, validated in manuscript
- **Input**: Pre-computed CSVs from S3 (65 MB, uploaded July 2024)
- **Processing**: Filtering and formatting only - NO recalculation
- **Output**: `data/processed/screen_results_baseline.duckdb` (178,826 rows)
- **Speed**: ~5 minutes total
- **Commands**: `just download-baseline && just run-baseline`
- **Code**: `haghighi_mito/data.py` (formatting only)

### Method 1: Regenerated - Notebook (Complete but Messy)
- **Status**: Complete pipeline, superseded by Method 2
- **Input**: Raw per-site profiles (2.7 GB)
- **Processing**: Calculate slopes + stats + filter + format
- **Output**: `data/processed/screen_results_notebook.duckdb`
- **Speed**: ~10 minutes per dataset
- **Commands**: `just download-raw && just run-notebook`
- **Code**: `notebooks/2.0-mh-virtual-screen.py` (1433 lines) + `data.py`
- **Note**: Exploratory code with dead branches (`if 0:`), but functional

### Method 2: Regenerated - Clean Module (Recommended for Development)
- **Status**: Complete pipeline, clean implementation
- **Input**: Raw per-site profiles (2.7 GB)
- **Processing**: Calculate slopes + stats + filter + format
- **Output**: `data/processed/screen_results_module.duckdb`
- **Speed**: ~10 minutes per dataset
- **Commands**: `just download-raw && just run-module`
- **Code**: `haghighi_mito/virtual_screen.py` (448 lines) + `diagnostics.py` + `data.py`
- **Note**: Professional refactor, identical output to Method 1, cleaner codebase

**See the Snakefile docstring for comprehensive pipeline documentation.**

## Critical Reproducibility Issue

**All regenerated methods (1 & 2) produce ~77% agreement with baseline due to lost July 2024 code.**

Key findings:
- Baseline generated with code that predates this repository (Sept 2025 creation)
- Input data confirmed identical (`Count_Cells_avg` matches 100%)
- 100% of large slope differences have different peak indices
- Root cause: Unknown peak detection algorithm used in July 2024
- Impact: Baseline trusted but opaque; regenerated methods transparent but inexact

**Match rates** (taorf, n=327):
- `t_target_pattern`: 37% within 10% (no peak detection, best metric)
- `slope`: 19% within 10% (peak detection differs)
- `t_orth`: 33% within 10% (validates processing)

**Recommendation**: Always use `just run-baseline` for validated/publication results. Use Methods 1/2 for development and future improvements.

See `docs/PROGRESS.md` for detailed investigation history.

## Dataset Configuration

Six datasets with specific metadata columns (defined in `haghighi_mito/config.py::DATASET_INFO`):

| Dataset | Key Column | Reference Set | Type |
|---------|-----------|---------------|------|
| LINCS | `Metadata_pert_name` | `Metadata_pert_id_dose` | Compounds |
| CDRP | `Metadata_pert_id` | `Metadata_Sample_Dose` | Compounds |
| JUMP-ORF | `Metadata_Symbol` | `Metadata_JCP2022` | Genetics |
| JUMP-CRISPR | `Metadata_Symbol` | `Metadata_JCP2022` | Genetics |
| JUMP-Compound | `Metadata_InChIKey` | `Metadata_JCP2022` | Compounds |
| TA-ORF | `Metadata_gene_name` | `Metadata_broad_sample` | Genetics |

## Virtual Screen Analysis

The pipeline identifies perturbations affecting mitochondrial radial distribution:

1. **Target Feature**: Primary metric is `slope` (radial distribution slope) or `d_slope` (effect size)
2. **Statistical Testing**: Per-site aggregation → per-plate tests vs. controls
3. **Filtering**:
   - Cell count filter: Remove bottom 10% by cell count
   - Target feature significance: Benjamini-Hochberg corrected p-values (FDR=0.05)
   - Orthogonal feature filter: Ensure perturbation doesn't affect unrelated features
4. **Output**: Ranked perturbation lists with statistical metrics

## Key Dependencies

Managed via Pixi (see `pyproject.toml`):

- **Pipeline**: snakemake, pixi, s5cmd (fast S3 downloads)
- **Data**: pandas, numpy, scipy, pyarrow (Parquet), duckdb
- **Visualization**: matplotlib, seaborn, scienceplots
- **Analysis**: sklearn, blitzgsea (enrichment)
- **Custom**: singlecell-morph package (morphological analysis)

## Important Notes

- **Data provenance**: Excel files in `curated_2024-08-11/` and `curated_2025-10-25/` are manually curated via Google Sheets, not direct pipeline outputs
- **Pipeline-generated files** are in `generated_from_*/` directories (gitignored, reproducible)
- **Always use** `pixi run python` instead of bare `python`
- **Never modify** `docs/PROGRESS.md` unless explicitly requested (see maintenance guidelines in that file)
- **Workflow principles**: Data flows one direction (raw → interim → processed), raw data is immutable
