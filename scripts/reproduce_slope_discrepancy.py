#!/usr/bin/env python3
"""
Slope Discrepancy Reproduction Script - TA-ORF Dataset

PURPOSE:
This script demonstrates the discrepancy between our current implementation and
the July 2024 baseline for the TA-ORF dataset.

KEY FINDING:
Imperfect reproduction of baseline slopes (r=0.849) and peak indices (r=0.516).

USAGE:
    # If running from this repository (pixi manages dependencies):
    pixi run python scripts/reproduce_slope_discrepancy.py

    # If running standalone (use uv to auto-manage dependencies):
    uv run --with pandas --with numpy --with scipy --with matplotlib \
        scripts/reproduce_slope_discrepancy.py

DEPENDENCIES:
    - pandas
    - numpy
    - scipy
    - matplotlib

WHAT IT DOES:
1. Loads data from data/external/mito_project/ (per-site profiles + metadata)
2. Calculates slopes using haghighi_mito/virtual_screen.py implementation
3. Compares with baseline from July 2024
4. Generates outputs:
   - Scatter plots showing discrepancy
   - CSV with detailed comparison
   - Summary statistics

FILES ACCESSED (for taorf dataset):
Input:
  - data/external/mito_project/workspace/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz
  - data/external/mito_project/workspace/per_site_aggregated_profiles_newpattern_2/taorf/2013_10_11_SIGMA2_Pilot_site_agg_profiles.csv.gz
  - data/external/mito_project/workspace/results/virtual_screen_baseline/taorf_results_pattern_aug_070624.csv
Output:
  - data/processed/virtual_screen_module/taorf_slope_discrepancy.png
  - data/processed/virtual_screen_module/taorf_slope_comparison.csv

EXPECTED RESULTS:
    Slope correlation: r = 0.849
    Peak index correlation: r = 0.516

CURRENT IMPLEMENTATION (haghighi_mito/virtual_screen.py):

Pipeline Steps:
1. Load data: Per-site aggregated profiles from S3
2. Filter plates: Keep only plates with control wells
3. Control subtraction: Subtract mean control profile per plate
4. Slope calculation:
   - Smooth 12-bin radial distribution (Savitzky-Golay filter, window_length=5, polyorder=3)
   - Find peaks: argmax and argmin of smoothed profile
   - Select last peak: max(argmax, argmin) excluding last 2 bins
   - Calculate slope: (mean(last 2 bins) - peak_value) / distance
5. Z-score normalization: Per-plate normalization of slopes
6. Aggregation: Median of all per-site slopes

Code Implementation:
See haghighi_mito/virtual_screen.py for:
- preprocess_metadata() - handles metadata loading and standardization
- load_dataset_data() - loads per-site profiles and merges with metadata
- calculate_metrics() - core slope calculation pipeline

QUESTIONS FOR ORIGINAL AUTHOR:

Given slope correlation r=0.849 and peak correlation r=0.516:

1. Is there additional preprocessing BEFORE control subtraction?
   - Code: haghighi_mito.virtual_screen.calculate_metrics()
   - Current: None (raw features used directly for control subtraction)
   - Note: Z-score normalization only applied to derived metrics (slope, last_peak_ind) AFTER calculation
   - Missing steps?

2. Is the smoothing window correct?
   - Code: haghighi_mito.vectorized_slope.find_end_slope2_vectorized()
   - Current: Savitzky-Golay filter, window_length=5, polyorder=3

3. Is the peak detection correct?
   - Code: haghighi_mito.vectorized_slope.find_end_slope2_vectorized()
   - Current: max(argmax(data), argmin(data)) excluding last 2 bins
   - Should this be different?

4. Is the slope endpoint calculation correct?
   - Code: haghighi_mito.vectorized_slope.find_end_slope2_vectorized()
   - Current: (mean(last 2 bins) - peak_value) / distance
   - Alternative: use only last bin?

5. Is the aggregation method correct?
   - Code: haghighi_mito.virtual_screen.calculate_metrics()
   - Current: median of ALL per-site slopes
   - Alternative: median of per-plate medians (two-stage)?

OUTPUTS:
After running, check:
    data/processed/virtual_screen_module/
    ├── taorf_slope_discrepancy.png      # Scatter plots (baseline vs current)
    └── taorf_slope_comparison.csv       # Detailed comparison for all 339 perturbations

BACKGROUND:
See docs/PROGRESS.md for full investigation history (Oct 25-27, 2025).

Key milestone: Oct 26 - Z-score normalization fix improved slope correlation from r=0 to r=0.849.

Remaining gap: Strong correlation but systematic differences suggest methodological
divergence, not just numerical precision issues.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from loguru import logger

# Install audit hook to trace file reads
files_read = []

def audit_hook(event, args):
    """Trace file open operations."""
    if event == 'open':
        filename = str(args[0])
        mode = args[1] if len(args) > 1 else ''
        # Only log data files
        if any(ext in filename for ext in ['.csv', '.parquet', '.xlsx', '.db']):
            files_read.append(filename)
            logger.info(f"Reading: {filename}")

sys.addaudithook(audit_hook)

# Add project root to path to import modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from haghighi_mito.config import DATASET_INFO
from haghighi_mito.virtual_screen import (
    calculate_metrics,
    load_dataset_data,
)


def main(
    dataset: str = typer.Argument(
        "taorf",
        help="Dataset to analyze (taorf, lincs, CDRP, jump_orf, jump_crispr, jump_compound)"
    )
):
    """Main execution."""
    logger.info("=" * 80)
    logger.info("MINIMAL SLOPE DISCREPANCY REPRODUCTION SCRIPT")
    logger.info("=" * 80)
    logger.info("This script runs our current slope calculation implementation")
    logger.info("and compares with the July 2024 baseline to identify discrepancies.")
    logger.info(f"Dataset: {dataset}")

    # Validate dataset
    if dataset not in DATASET_INFO:
        logger.error(f"Unknown dataset: {dataset}")
        logger.error(f"Valid datasets: {', '.join(DATASET_INFO.keys())}")
        sys.exit(1)

    # Check that data exists
    mito_project_root = Path("data/external/mito_project")
    data_dir = mito_project_root / "workspace"
    baseline_file = data_dir / f"results/virtual_screen_baseline/{dataset}_results_pattern_aug_070624.csv"

    if not baseline_file.exists():
        logger.error("Baseline data not found!")
        logger.error(f"Expected: {baseline_file}")
        logger.error("Please run first: just generate-baseline-all")
        sys.exit(1)

    logger.info("=" * 80)
    logger.info("STEP 1: LOAD DATA")
    logger.info("=" * 80)

    # Load data using module functions (which internally calls preprocess_metadata)
    per_site_df, annot = load_dataset_data(dataset)

    # Get dataset-specific configuration (early, for logging)
    pert_col = DATASET_INFO[dataset]["pert_col"]

    logger.info(f"Loaded data:")
    logger.info(f"  - Per-site profiles: {len(per_site_df)} rows")
    logger.info(f"  - Metadata: {len(annot)} perturbations")
    logger.info(f"  - Unique perturbations: {per_site_df[pert_col].nunique()}")

    logger.info("=" * 80)
    logger.info("STEP 2: CALCULATE SLOPES (Our Implementation)")
    logger.info("=" * 80)
    logger.info("Pipeline steps:")
    logger.info("  1. Filter to plates with controls")
    logger.info("  2. Subtract control mean per plate from radial distribution")
    logger.info("  3. Calculate slopes on control-subtracted profiles (vectorized)")
    logger.info("  4. Z-score normalize slopes per plate")
    logger.info("  5. Aggregate via median across plates per perturbation")

    # Calculate slopes using module function
    results, per_site_df_with_slopes = calculate_metrics(per_site_df, annot, dataset)

    logger.info(f"Calculated slopes for {len(results)} perturbations")

    # Get dataset-specific metadata columns
    meta_cols = DATASET_INFO[dataset]["meta_cols"]

    # Keep only the columns we need for comparison
    comparison_cols = meta_cols + ["Count_Cells_avg", "slope", "last_peak_ind"]
    current = results[comparison_cols].copy()

    logger.info("=" * 80)
    logger.info("STEP 3: COMPARE WITH BASELINE")
    logger.info("=" * 80)

    # Load baseline
    baseline = pd.read_csv(baseline_file)
    logger.info(f"Loaded baseline: {len(baseline)} perturbations")

    # Check that required columns exist in baseline
    missing_cols = [col for col in comparison_cols if col not in baseline.columns]
    if missing_cols:
        logger.error(f"Baseline file missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(baseline.columns)}")
        sys.exit(1)

    # Drop duplicates in pert_col to avoid Cartesian product during merge
    # (e.g., DMSO controls with multiple metadata values)
    baseline_dedup = baseline[comparison_cols].drop_duplicates(
        subset=pert_col, keep="first"
    )
    current_dedup = current[comparison_cols].drop_duplicates(
        subset=pert_col, keep="first"
    )

    n_baseline_dropped = len(baseline) - len(baseline_dedup)
    n_current_dropped = len(current) - len(current_dedup)

    if n_baseline_dropped > 0:
        logger.warning(f"Dropped {n_baseline_dropped} duplicate {pert_col} from baseline ({len(baseline_dedup)} remaining)")
    if n_current_dropped > 0:
        logger.warning(f"Dropped {n_current_dropped} duplicate {pert_col} from current ({len(current_dedup)} remaining)")

    # Merge on the perturbation column
    merge_cols = [pert_col, "Count_Cells_avg", "slope", "last_peak_ind"]
    comparison = baseline_dedup.merge(
        current_dedup[merge_cols],
        on=pert_col,
        suffixes=("_baseline", "_current"),
        how="inner"
    )

    logger.info(f"Matched: {len(comparison)} perturbations")

    # Calculate metrics
    corr_cells = comparison[["Count_Cells_avg_baseline", "Count_Cells_avg_current"]].corr().iloc[0, 1]
    corr_slope = comparison[["slope_baseline", "slope_current"]].corr().iloc[0, 1]
    corr_peak = comparison[["last_peak_ind_baseline", "last_peak_ind_current"]].corr().iloc[0, 1]

    comparison["slope_pct_diff"] = (
        100 * abs(comparison["slope_current"] - comparison["slope_baseline"]) /
        (abs(comparison["slope_baseline"]) + 1e-10)
    )

    comparison["peak_diff"] = abs(
        comparison["last_peak_ind_current"] - comparison["last_peak_ind_baseline"]
    )

    logger.info("=" * 80)
    logger.info("RESULTS SUMMARY")
    logger.info("=" * 80)
    logger.info(f"  Cell count correlation: r = {corr_cells:.3f}")
    logger.info(f"  Slope correlation: r = {corr_slope:.3f}")
    logger.info(f"  Peak index correlation: r = {corr_peak:.3f}")
    logger.info("Detailed comparison saved to CSV for inspection.")

    # Generate plots
    logger.info("=" * 80)
    logger.info("GENERATING PLOTS")
    logger.info("=" * 80)

    _fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Cell count comparison
    ax = axes[0]
    ax.scatter(comparison["Count_Cells_avg_baseline"], comparison["Count_Cells_avg_current"],
               alpha=0.5, s=30)

    all_cells = pd.concat([comparison["Count_Cells_avg_baseline"], comparison["Count_Cells_avg_current"]])
    cells_range = [all_cells.min(), all_cells.max()]
    ax.plot(cells_range, cells_range, 'r--', alpha=0.5, label="y=x")

    ax.set_xlabel("Baseline cell count")
    ax.set_ylabel("Current implementation cell count")
    ax.set_title(f"Cell Count Comparison (r={corr_cells:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Slope comparison
    ax = axes[1]
    ax.scatter(comparison["slope_baseline"], comparison["slope_current"],
               alpha=0.5, s=30)

    # Add diagonal line (perfect agreement)
    all_slopes = pd.concat([comparison["slope_baseline"], comparison["slope_current"]])
    slope_range = [all_slopes.min(), all_slopes.max()]
    ax.plot(slope_range, slope_range, 'r--', alpha=0.5, label="y=x")

    ax.set_xlabel("Baseline slope (July 2024)")
    ax.set_ylabel("Current implementation slope")
    ax.set_title(f"Slope Comparison (r={corr_slope:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Peak comparison
    ax = axes[2]
    ax.scatter(comparison["last_peak_ind_baseline"], comparison["last_peak_ind_current"],
               alpha=0.5, s=30)

    # Calculate dynamic range for z-scored peak indices (not raw bins!)
    all_peaks = pd.concat([comparison["last_peak_ind_baseline"], comparison["last_peak_ind_current"]])
    peak_range = [all_peaks.min(), all_peaks.max()]
    ax.plot(peak_range, peak_range, 'r--', alpha=0.5, label="y=x")

    ax.set_xlabel("Baseline peak index")
    ax.set_ylabel("Current implementation peak index")
    ax.set_title(f"Peak Index Comparison (r={corr_peak:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_dir = Path("data/processed/virtual_screen_module")
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_file = output_dir / f"{dataset}_slope_discrepancy.png"
    plt.savefig(plot_file, dpi=150)
    logger.success(f"Saved plot: {plot_file}")

    # Save comparison CSV
    csv_file = output_dir / f"{dataset}_slope_comparison.csv"
    comparison.to_csv(csv_file, index=False)
    logger.success(f"Saved detailed comparison: {csv_file}")

    # Summary
    logger.info("=" * 80)
    logger.info("NOTES")
    logger.info("=" * 80)
    logger.info("For implementation details and questions for the original author,")
    logger.info("see the docstring at the top of this file:")
    logger.info("  python -c \"import scripts.reproduce_slope_discrepancy; help(scripts.reproduce_slope_discrepancy)\"")
    logger.info("Or just read the docstring in the source code.")

    # Print file trace summary
    logger.info("=" * 80)
    logger.info("FILE TRACE SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total files accessed: {len(files_read)}")
    logger.info(f"Unique files: {len(set(files_read))}")
    for f in sorted(set(files_read)):
        logger.info(f"  - {f}")

    logger.info("=" * 80)
    logger.info("SCRIPT COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Outputs saved to: {output_dir}/")


if __name__ == "__main__":
    typer.run(main)
