# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis pipelines and results for a mitochondrial morphology study using Cell Painting assays. The project performs virtual screens across multiple perturbation datasets (compounds, ORFs, CRISPR) to identify genetic and chemical modulators of mitochondrial radial distribution patterns in fibroblasts.

## Repository Structure

This repository follows the [Carpenter-Singh lab workflow conventions](protocols/workflows.md) based on Cookiecutter Data Science principles:

```
2025_Haghighi_Mito/
├── data/                  # Project data following lab conventions
│   ├── external/          # Third-party reference datasets (KEGG, WikiPathways)
│   ├── interim/           # Intermediate processing outputs
│   └── processed/         # Final analysis outputs
│       ├── figures/       # Publication figures (PDFs)
│       └── tables/        # Screen results (XLSX)
├── notebooks/             # Analysis notebooks (numbered by phase)
├── pipelines/             # CellProfiler pipelines (.cppipe files)
├── protocols/             # Experimental protocols and workflow documentation
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

Notebooks follow Carpenter-Singh lab convention: `<phase>.<sequence>-<initials>-<description>.ipynb`

- Phase 1 = Data cleaning/feature engineering
- Phase 2 = Analysis/screening

- **`1.0-mh-feat-importance.ipynb`**: Feature importance analysis for mitochondrial morphology measurements.

- **`2.0-mh-virtual-screen.ipynb`**: Primary pipeline for screening phenotype strength of mitochondrial radial distribution features across datasets (LINCS, CDRP, JUMP-ORF, JUMP-CRISPR, JUMP-Compound, TA-ORF). Performs statistical testing, filtering, and generates ranked perturbation lists.

- **`2.1-mh-set-enrichment-analysis.ipynb`**: Gene set enrichment analysis (GSEA) using blitzgsea package. Analyzes enrichment in GO terms, pathways (KEGG, WikiPathways), disease associations (OMIM, HPO), and mechanism-of-action (MOA) categories.

- **`2.2-mh-check-vs-lists.ipynb`**: Validation and comparison of virtual screen hit lists.

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

- **`data/external/`**: Third-party reference datasets
  - `KEGG_2021_Human_table.txt`: KEGG pathway annotations for enrichment analysis
  - `WikiPathways_2024_Human_table.txt`: WikiPathways annotations for enrichment analysis

- **`data/interim/`**: Intermediate processing outputs
  - Currently empty (large profile datasets remain at remote locations)

- **`data/processed/`**: Final analysis outputs from notebooks
  - **`data/processed/tables/`**: Virtual screen results for each dataset
    - `CDRP_screen_results.xlsx`: CDRP compound screen results
    - `jump_compound_screen_results.xlsx`: JUMP compound screen results
    - `jump_crispr_screen_results.xlsx`: JUMP CRISPR screen results
    - `jump_orf_screen_results.xlsx`: JUMP ORF screen results
    - `lincs_screen_results.xlsx`: LINCS compound screen results
    - `taorf_screen_results.xlsx`: TA-ORF screen results
  - **`data/processed/figures/`**: Publication figures
    - `Figure3.pdf`, `Figure4b.pdf`, `Figure4c.pdf`: Main figures
    - `SuppFigure2.pdf`: Supplementary figure
    - `dendrogram_mito.pdf`: Hierarchical clustering dendrogram

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

### Jupyter Notebooks

All analysis is performed in Jupyter notebooks (`.ipynb` files). These notebooks:

- Connect to remote data sources (S3 buckets at paths like `/home/jupyter-mhaghigh@broadinst-ee45a/bucket/`)
- Use custom single-cell analysis package: `singlecell` (imported from `SingleCell_Morphological_Analysis/`)
- Read external reference data from `data/external/`
- Generate figures saved to `data/processed/figures/`
- Generate screen results saved to `data/processed/tables/`

### CellProfiler Pipelines

- Pipelines use `.cppipe` format (CellProfiler pipeline files)
- Batch processing versions include `_batch` suffix
- CP3-compatible versions include `_cp3` suffix
- To modify pipelines, open in CellProfiler application

### Key Python Dependencies

- pandas, numpy, scipy: Data analysis
- matplotlib, seaborn: Visualization
- sklearn: Preprocessing and normalization
- blitzgsea: Gene set enrichment analysis
- Custom `singlecell` package for morphological analysis

### Important File Paths

Remote data paths referenced in notebooks (may not be accessible locally):

- Profiles: `~/gallery/cpg0016-jump/`, `~/gallery/cpg0012-wawer-bioactivecompoundprofiling/`
- Project workspace: `~/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/`

## Analysis Workflow

1. **CellProfiler Processing**: Run `.cppipe` pipelines on microscopy images to extract morphological features
2. **Feature Analysis**: Use `notebooks/1.0-mh-feat-importance.ipynb` to analyze feature importance
3. **Virtual Screening**: Use `notebooks/2.0-mh-virtual-screen.ipynb` to identify hits across datasets
4. **Enrichment Analysis**: Use `notebooks/2.1-mh-set-enrichment-analysis.ipynb` to find enriched pathways/gene sets
5. **Validation**: Use `notebooks/2.2-mh-check-vs-lists.ipynb` to compare and validate results

## Critical Analysis Parameters

- **BH-corrected critical values** (FDR=0.05): Dataset-specific, stored in `orth_bh_corrected_critical_dict` and `target_bh_corrected_critical_dict`
- **Minimum gene set size**: 4-10 depending on analysis
- **Cell count quantile filter**: 0.1 (bottom 10%)
- **Target feature**: `slope` or `d_slope` (effect size of radial distribution slope)

## Important Instructions for Claude Code

- Do what has been asked; nothing more, nothing less
- NEVER create files unless absolutely necessary for achieving your goal
- ALWAYS prefer editing an existing file to creating a new one
- NEVER proactively create documentation files (*.md) or README files unless explicitly requested
- NEVER use emojis in any output, documentation, code comments, or commit messages unless explicitly requested
- Always use `pixi run python` instead of just `python`