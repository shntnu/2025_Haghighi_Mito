# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis pipelines and results for a mitochondrial morphology study using Cell Painting assays. The project performs virtual screens across multiple perturbation datasets (compounds, ORFs, CRISPR) to identify genetic and chemical modulators of mitochondrial radial distribution patterns in fibroblasts.

Supporting repository for:
> Haghighi, M. et al. Identifying and targeting abnormal mitochondrial localization associated with psychoses. *bioRxiv* 2025.10.08.676630 (2025) doi:[10.1101/2025.10.08.676630](https://doi.org/10.1101/2025.10.08.676630).

## Repository Structure

This repository follows the [Carpenter-Singh lab workflow conventions](protocols/workflows.md):

```text
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

```bash
just --list                          # View all available commands
```

**Quick workflows:**

```bash
# Baseline (validated, fast ~5 min)
just generate-baseline-all           # Download from S3 → Excel → DuckDB

# Regenerate from raw data (~10 min/dataset)
just generate-module-all             # Full pipeline: profiles → CSV → Excel → DuckDB
just generate-module-csv-for taorf   # Single dataset CSV only

# Diagnostics (~1 sec)
just diagnose-for taorf              # Compare with baseline + plots
just diagnose-all                    # All datasets
```

**Python CLI (for advanced control):**

```bash
pixi run haghighi-mito --help
pixi run haghighi-mito virtual-screen --dataset taorf
pixi run haghighi-mito compare-baseline --dataset taorf
```

**Development:**

```bash
pixi run python / jupyter lab / pytest / ruff check .
```

## Pipeline Architecture

Three methods for generating virtual screen results:

| Method | Input | Processing | Speed | Command | Use Case |
|--------|-------|------------|-------|---------|----------|
| **0: Baseline** | Pre-computed CSVs (65 MB) | Filtering only | ~5 min | `just generate-baseline-all` | **Production** (validated) |
| **1: Notebook** | Raw profiles (2.7 GB) | Full recalculation | ~10 min/dataset | `just generate-notebook-all` | Reference (original code) |
| **2: Module** | Raw profiles (2.7 GB) | Full recalculation | ~10 min/dataset | `just generate-module-all` | **Development** (clean code) |

- **Method 0** outputs: `screen_results_baseline.duckdb` (178,826 rows), uses `haghighi_mito/data.py`
- **Method 1** outputs: `screen_results_notebook.duckdb`, uses `notebooks/2.0-mh-virtual-screen.py` (1433 lines)
- **Method 2** outputs: `screen_results_module.duckdb`, uses `haghighi_mito/virtual_screen.py` (448 lines) + `diagnostics.py`

See Snakefile docstring for full pipeline documentation.

## Reproducibility Status

⚠️ **Regenerated methods (1 & 2) show incomplete agreement with validated baseline.**

- **Root cause**: Baseline generated with pre-repository code (Sept 2025 repo creation)
- **Input data**: Confirmed identical (100% match on `Count_Cells_avg`)
- **Agreement**: Method 2 closer to baseline than Method 1, but discrepancies remain unexplained
- **Recommendation**: Use Method 0 for publication; Method 2 for development

See `docs/PROGRESS.md` for detailed diagnostics and investigation history.

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
