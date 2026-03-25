# 2025_Haghighi_Mito - Reproducible Workflow Fork

This fork investigated whether the virtual screening analysis from the paper could be cleanly reproduced and implemented as a production-quality workflow.

Supporting paper: Haghighi, M. et al. Identifying and targeting abnormal mitochondrial localization associated with psychoses. *bioRxiv* 2025.10.08.676630 (2025) doi:[10.1101/2025.10.08.676630](https://doi.org/10.1101/2025.10.08.676630).

## What we found

The analysis is fully reproducible. We identified two algorithmic fixes required for exact baseline reproduction:

1. Pre-standardize radial features per plate before control subtraction
2. Use `nanpercentile(interpolation="nearest")` for median plate selection

With these fixes, the clean module implementation (`haghighi_mito/virtual_screen.py`) achieved perfect correlation (r=1.000) with the July 2024 baseline for taorf and lincs datasets.

See `docs/PROGRESS.md` for detailed investigation history.

## Status

The reproduction work is complete. The upstream paper repository at [carpenter-singh-lab/2025_Haghighi_Mito](https://github.com/carpenter-singh-lab/2025_Haghighi_Mito) continues active development. Upstream changes are tracked via a git worktree at `.worktrees/upstream/`.
