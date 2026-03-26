# 2025_Haghighi_Mito

Reproducible analysis pipelines for mitochondrial morphology in patient fibroblasts using Cell Painting assays. Performs virtual screens across multiple perturbation datasets (compounds, ORFs, CRISPR) to identify genetic and chemical modulators of mitochondrial radial distribution patterns.

Supporting paper: Haghighi, M. et al. Identifying and targeting abnormal mitochondrial localization associated with psychoses. *bioRxiv* 2025.10.08.676630 (2025) doi:[10.1101/2025.10.08.676630](https://doi.org/10.1101/2025.10.08.676630).

## What this repo contains

**Virtual screen pipeline** (`haghighi_mito/virtual_screen.py`): Tested, vectorized implementation that reproduces the July 2024 baseline results (r=1.000 for all metrics). Processes 6 datasets covering ~150K perturbations.

**Patient phenotype analysis** (`haghighi_mito/patient_analysis.py` + `notebooks/1.0-*.py`, `3.0-*.py`): Single-cell QC, feature standardization, and MITO-SLOPE computation from 185 patient fibroblast SQLite databases. Generates Figure 3, 4b, 4c, SuppFigure 2, and dendrogram.

**Supplemental figure generation** (`notebooks/3.1-*.py`, `3.2-*.py`): 16 intensity/skeleton feature boxplots and 24 radial distribution polar plots, all reproducible from a single pre-aggregated CSV.

**Enrichment analysis** (`notebooks/2.1-*.py`): GSEA using blitzgsea across KEGG, WikiPathways, and custom gene/drug sets.

All manuscript figures except Figure 2 (manual microscopy composite) are reproducible from this repository.

## Quick start

```bash
# Virtual screen (validated, ~5 min)
snakemake all_baseline -c4 -p

# Regenerate from raw data (~10 min/dataset)
snakemake all_module -c4 -p

# Patient phenotype figures
snakemake download_all_patient_data -c1 -p   # Download ~2.9 GB from S3
pixi run python notebooks/1.0-mh-feat-importance.py  # Figure 3, 4b, SuppFig 2, dendrogram
pixi run python notebooks/3.0-mh-slope-analysis.py   # Figure 4c
pixi run python notebooks/3.1-ew-supplemental-mito-features.py  # 16 supplemental plots
pixi run python notebooks/3.2-ew-radial-distribution-plots.py   # 24 polar plots
```

See `CLAUDE.md` for full command reference and `docs/PROGRESS.md` for investigation history.

## Relationship to upstream

This repo originated as a fork of [carpenter-singh-lab/2025_Haghighi_Mito](https://github.com/carpenter-singh-lab/2025_Haghighi_Mito). All upstream notebook logic has been captured here — either distilled into tested Python modules or converted to runnable jupytext notebooks with config-based paths. Upstream changes are tracked via a git worktree at `.worktrees/upstream/`.
