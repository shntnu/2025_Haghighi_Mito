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
│   ├── virtual_screen.py  # Clean module: per-site profiles → virtual screen CSVs (~656 lines)
│   ├── diagnostics.py     # Baseline comparison plots and metrics
│   ├── vectorized_slope.py # Optimized slope calculation (~200x speedup)
│   ├── vectorized_stats.py # Batch statistical testing
│   └── tests/             # Unit tests (test_vectorized_slope.py, test_vectorized_stats.py)
├── notebooks/             # Jupytext .py notebooks (not .ipynb — use jupytext or marimo to run)
│   ├── 1.0-mh-feat-importance.py           # Feature importance (patient fibroblasts)
│   ├── 2.0-mh-virtual-screen.py            # Virtual screening (1438 lines)
│   ├── 2.0-mh-virtual-screen-original.py   # Original baseline version (archived)
│   ├── 2.0-mh-virtual-screen-minimal.py    # Minimal version for testing
│   ├── 2.1-mh-set-enrichment-analysis.py   # GSEA (blitzgsea)
│   ├── 2.2-mh-check-vs-lists.py            # Validation and filtering
│   └── explore.py                           # Exploratory analysis
├── data/
│   ├── external/          # Downloaded data from S3 (~5.1 GB, gitignored)
│   ├── interim/           # Intermediate Parquet files (gitignored)
│   └── processed/         # Final outputs (DuckDB, Excel, figures)
├── scripts/               # Utility scripts
│   ├── reproduce_slope_discrepancy.py  # Diagnostic: reproduce slopes from raw profiles
│   ├── restore_intelligent.py          # S3 Glacier data restoration
│   ├── generate-dag.sh                 # Pipeline DAG visualization
│   └── check_restore_status.sh         # S3 restore status checker
├── docs/                  # Project documentation
│   ├── PROGRESS.md        # Investigation history and milestones (DO NOT MODIFY)
│   ├── DATA_FLOW.md       # Complete data flow architecture
│   ├── DATA_DOWNLOAD.md   # Data requirements (~5.1 GB breakdown)
│   ├── MANUSCRIPT.md      # Manuscript preparation notes
│   └── pipeline/          # DAG visualizations (PNG)
├── pipelines/             # CellProfiler pipeline XMLs (15+ files, image analysis)
├── protocols/             # Experimental protocols and workflow docs
├── tests/                 # Standalone tests (scipy_version_test)
├── Snakefile              # Pipeline automation (all commands)
├── pyproject.toml         # Pixi configuration and dependencies
└── flake.nix              # NixOS environment (+ .envrc for direnv auto-activation)
```

## Common Commands

```bash
snakemake --list-target-rules        # View all available targets
```

**Quick workflows:**

```bash
# Baseline (validated, fast ~5 min)
snakemake all_baseline -c4 -p        # Download from S3 → Excel → DuckDB

# Regenerate from raw data (~10 min/dataset)
snakemake all_module -c4 -p          # Full pipeline: profiles → CSV → Excel → DuckDB

# Single dataset CSV only
snakemake data/processed/virtual_screen_module/taorf_results_pattern_aug_070624.csv -c1 -p

# Diagnostics (~1 sec)
snakemake data/processed/virtual_screen_module/taorf_comparison_metrics.png -c1 -p  # Single
snakemake all_module_diagnose -c4 -p                                                 # All

# Full analysis (diagnostics + reproduction scripts)
snakemake all_analysis -c4 -p
```

**Utility commands:**

```bash
snakemake -n -p <target>             # Dry run (preview what will run)
snakemake --summary                  # Show pipeline status
snakemake generate_dag_visualizations  # Generate DAG PNGs
snakemake validate_database_pair --config db1=baseline db2=module  # Compare databases
```

**Python CLI (for advanced control):**

```bash
pixi run haghighi-mito --help
pixi run haghighi-mito virtual-screen --dataset taorf
pixi run haghighi-mito compare-baseline --dataset taorf
```

**Testing:**

```bash
pixi run pytest                      # Run all tests (haghighi_mito/tests/)
pixi run pytest -x -v                # Stop on first failure, verbose
```

Tests cover vectorized slope calculation (9 tests) and vectorized statistics (16 tests).

**Development:**

```bash
pixi run python / jupyter lab / pytest / ruff check .
```

## Pipeline Architecture

Three methods for generating virtual screen results:

| Method | Input | Processing | Speed | Command | Use Case |
|--------|-------|------------|-------|---------|----------|
| **0: Baseline** | Pre-computed CSVs (65 MB) | Filtering only | ~5 min | `snakemake all_baseline` | **Production** (validated) |
| **1: Notebook** | Raw profiles (2.7 GB) | Full recalculation | ~10 min/dataset | `snakemake all_notebook` | Reference (original code) |
| **2: Module** | Raw profiles (2.7 GB) | Full recalculation | ~10 min/dataset | `snakemake all_module` | **Development** (clean code) |

- **Method 0** outputs: `screen_results_baseline.duckdb` (178,826 rows), uses `haghighi_mito/data.py`
- **Method 1** outputs: `screen_results_notebook.duckdb`, uses `notebooks/2.0-mh-virtual-screen.py` (1438 lines)
- **Method 2** outputs: `screen_results_module.duckdb`, uses `haghighi_mito/virtual_screen.py` (~656 lines) + `diagnostics.py`

See Snakefile docstring for full pipeline documentation.

## Reproducibility Status

✅ **Module perfectly reproduces July 2024 baseline** (r=1.000 for all metrics).

Two algorithmic fixes resolved prior discrepancies:
1. Pre-standardize radial features per plate before control subtraction
2. Use `nanpercentile` for median plate selection (matching upstream method)

- **Validation**: taorf (324 perturbations) and lincs (9,395 perturbations) both show 99.7-100% agreement
- **Recommendation**: Use Method 2 (module) for both development and production

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

## Code Style

Ruff with permissive settings (see `pyproject.toml`):
- **Line length**: 300 (intentionally wide — scientific code with long variable names)
- **Target**: Python 3.12
- Notebooks have relaxed rules (E402, F401, E703 etc. are ignored)

## Git Workflow

- **origin**: `shntnu/2025_Haghighi_Mito` (personal fork)
- **upstream**: `carpenter-singh-lab/2025_Haghighi_Mito` (lab repo)
- **Upstream worktree**: `.worktrees/upstream/` tracks `upstream/main`. Update with `git fetch upstream` then `cd .worktrees/upstream && git pull`

## Important Notes

- **Data provenance**: Excel files in `data/external/tables/curated_2024-08-11/` and `data/external/tables/curated_2025-10-25/` are manually curated via Google Sheets (originate from upstream `carpenter-singh-lab/2025_Haghighi_Mito`), not direct pipeline outputs
- **Pipeline-generated files** are in `generated_from_*/` directories (gitignored, reproducible)
- **Always use** `pixi run python` instead of bare `python`
- **Never modify** `docs/PROGRESS.md` unless explicitly requested (see maintenance guidelines in that file)
- **Workflow principles**: Data flows one direction (raw → interim → processed), raw data is immutable
- **Notebooks are Jupytext `.py` files**, not `.ipynb`. They have YAML frontmatter for Jupyter metadata. Edit as Python files; run with `pixi run jupyter lab` (jupytext syncs automatically) or `pixi run marimo edit`
