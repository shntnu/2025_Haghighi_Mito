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

EXPECTED RESULTS:
    Slope correlation: r = 0.849
    Peak index correlation: r = 0.516

CURRENT IMPLEMENTATION (haghighi_mito/virtual_screen.py):

Pipeline Steps:
1. Load data: Per-site aggregated profiles from S3
2. Filter plates: Keep only plates with control wells
3. Control subtraction: Subtract mean control profile per plate
4. Slope calculation:
   - Smooth 12-bin radial distribution (moving average, window=3)
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
   - Current: z-score per plate on radial + orthogonal features
   - Missing steps?

2. Is the smoothing window correct?
   - Current: moving average, window size = 3

3. Is the peak detection correct?
   - Current: max(argmax(data), argmin(data)) excluding last 2 bins
   - Should this be different?

4. Is the slope endpoint calculation correct?
   - Current: (mean(last 2 bins) - peak_value) / distance
   - Alternative: use only last bin?

5. Is the aggregation method correct?
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

# Add project root to path to import modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from haghighi_mito.virtual_screen import (
    calculate_metrics,
    load_dataset_data,
)


def main():
    """Main execution."""
    print("\n" + "=" * 80)
    print("MINIMAL SLOPE DISCREPANCY REPRODUCTION SCRIPT")
    print("=" * 80)
    print("\nThis script runs our current slope calculation implementation")
    print("and compares with the July 2024 baseline to identify discrepancies.\n")

    dataset = "taorf"

    # Check that data exists
    mito_project_root = Path("data/external/mito_project")
    data_dir = mito_project_root / "workspace"
    baseline_file = data_dir / "results/virtual_screen_baseline/taorf_results_pattern_aug_070624.csv"

    if not baseline_file.exists():
        print("ERROR: Baseline data not found!")
        print(f"Expected: {baseline_file}")
        print("\nPlease run first: just generate-baseline-all")
        sys.exit(1)

    print("=" * 80)
    print("STEP 1: LOAD DATA")
    print("=" * 80)

    # Load data using module functions (which internally calls preprocess_metadata)
    per_site_df, annot = load_dataset_data(dataset)

    print(f"\nLoaded data:")
    print(f"  - Per-site profiles: {len(per_site_df)} rows")
    print(f"  - Metadata: {len(annot)} perturbations")
    print(f"  - Unique perturbations: {per_site_df['Metadata_broad_sample'].nunique()}")

    print("\n" + "=" * 80)
    print("STEP 2: CALCULATE SLOPES (Our Implementation)")
    print("=" * 80)
    print("\nPipeline steps:")
    print("  1. Filter to plates with controls")
    print("  2. Subtract control mean per plate from radial distribution")
    print("  3. Calculate slopes on control-subtracted profiles (vectorized)")
    print("  4. Z-score normalize slopes per plate")
    print("  5. Aggregate via median across plates per perturbation\n")

    # Calculate slopes using module function
    results, per_site_df_with_slopes = calculate_metrics(per_site_df, annot, dataset)

    print(f"\nCalculated slopes for {len(results)} perturbations")

    # Keep only the columns we need for comparison
    current = results[["Metadata_broad_sample", "Metadata_gene_name", "slope", "last_peak_ind"]].copy()

    print("\n" + "=" * 80)
    print("STEP 3: COMPARE WITH BASELINE")
    print("=" * 80)

    # Load baseline
    baseline = pd.read_csv(baseline_file)
    print(f"\nLoaded baseline: {len(baseline)} perturbations")

    # Merge
    comparison = baseline[["Metadata_broad_sample", "Metadata_gene_name", "slope", "last_peak_ind"]].merge(
        current[["Metadata_broad_sample", "slope", "last_peak_ind"]],
        on="Metadata_broad_sample",
        suffixes=("_baseline", "_current"),
        how="inner"
    )

    print(f"Matched: {len(comparison)} perturbations")

    # Calculate metrics
    corr_slope = comparison[["slope_baseline", "slope_current"]].corr().iloc[0, 1]
    corr_peak = comparison[["last_peak_ind_baseline", "last_peak_ind_current"]].corr().iloc[0, 1]

    comparison["slope_pct_diff"] = (
        100 * abs(comparison["slope_current"] - comparison["slope_baseline"]) /
        (abs(comparison["slope_baseline"]) + 1e-10)
    )

    comparison["peak_diff"] = abs(
        comparison["last_peak_ind_current"] - comparison["last_peak_ind_baseline"]
    )

    print(f"\n{'=' * 80}")
    print("RESULTS SUMMARY")
    print(f"{'=' * 80}")
    print(f"  Slope correlation: r = {corr_slope:.3f}")
    print(f"  Peak index correlation: r = {corr_peak:.3f}")
    print(f"\nDetailed comparison saved to CSV for inspection.")

    # Generate plots
    print(f"\n{'=' * 80}")
    print("GENERATING PLOTS")
    print(f"{'=' * 80}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Slope comparison
    ax = axes[0]
    ax.scatter(comparison["slope_baseline"], comparison["slope_current"],
               alpha=0.5, s=30)

    # Add diagonal line (perfect agreement)
    all_slopes = pd.concat([comparison["slope_baseline"], comparison["slope_current"]])
    slope_range = [all_slopes.min(), all_slopes.max()]
    ax.plot(slope_range, slope_range, 'r--', alpha=0.5, label="Perfect agreement")

    ax.set_xlabel("Baseline slope (July 2024)")
    ax.set_ylabel("Current implementation slope")
    ax.set_title(f"Slope Comparison (r={corr_slope:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Peak comparison
    ax = axes[1]
    ax.scatter(comparison["last_peak_ind_baseline"], comparison["last_peak_ind_current"],
               alpha=0.5, s=30)

    peak_range = [0, 12]
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
    print(f"\n✓ Saved plot: {plot_file}")

    # Save comparison CSV
    csv_file = output_dir / f"{dataset}_slope_comparison.csv"
    comparison.to_csv(csv_file, index=False)
    print(f"✓ Saved detailed comparison: {csv_file}")

    # Summary
    print("\n" + "=" * 80)
    print("NOTES")
    print("=" * 80)
    print("For implementation details and questions for the original author,")
    print("see the docstring at the top of this file:")
    print("  python -c \"import scripts.reproduce_slope_discrepancy; help(scripts.reproduce_slope_discrepancy)\"")
    print("\nOr just read the docstring in the source code.")

    print("=" * 80)
    print("SCRIPT COMPLETE")
    print("=" * 80)
    print(f"\nOutputs saved to: {output_dir}/")


if __name__ == "__main__":
    main()
