"""Compare regenerated (module) results with validated baseline (July 2024).

Functions:
    compare_with_baseline_csv(dataset) - Main entry: load CSVs, compute diffs, generate plots
    plot_baseline_comparison(dataset) - Internal: create 2x2 scatter plots (auto-called)
    compare_with_baseline(results, dataset) - Internal: comparison logic

Workflow:
    just diagnose-for taorf  → comparison CSV + diagnostic PNG (1 command!)

Files:
    Input:  virtual_screen_module/{dataset}_results_pattern_aug_070624.csv (module)
            virtual_screen_baseline/{dataset}_results_pattern_aug_070624.csv (S3)
    Output: virtual_screen_module/{dataset}_baseline_comparison.csv (comparison)
            figures/diagnostics/{dataset}_comparison_metrics.png (plots)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import linregress, pearsonr, spearmanr

from haghighi_mito.config import DATASET_INFO, EXTERNAL_DATA_DIR, PROCESSED_DATA_DIR


def compare_with_baseline(results, dataset: str):
    """Load baseline CSV, merge with results, calculate diffs (internal helper).

    Returns DataFrame with columns: *_new, *_baseline, *_diff, *_pct_diff
    Prints summary statistics to logger.
    """
    baseline_path = EXTERNAL_DATA_DIR / "mito_project/workspace/results/virtual_screen_baseline" / f"{dataset}_results_pattern_aug_070624.csv"

    logger.info(f"Loading baseline from {baseline_path}")
    baseline = pd.read_csv(baseline_path)

    pert_col = DATASET_INFO[dataset]["pert_col"]
    baseline_cols = [pert_col, "Count_Cells_avg", "last_peak_ind", "slope"]

    # Add t-values if they exist
    if "t_target_pattern" in results.columns:
        baseline_cols.extend(["t_target_pattern", "t_orth", "t_slope", "d_slope"])

    comparison = pd.merge(results, baseline[baseline_cols], on=pert_col, suffixes=("_new", "_baseline"))

    logger.info(f"Matched {len(comparison)} perturbations between new and baseline")

    # Calculate differences: absolute (_diff) and percentage (_pct_diff)
    comparison["Count_Cells_diff"] = comparison["Count_Cells_avg_new"] - comparison["Count_Cells_avg_baseline"]
    comparison["Count_Cells_pct_diff"] = 100 * comparison["Count_Cells_diff"] / comparison["Count_Cells_avg_baseline"]

    comparison["last_peak_ind_diff"] = comparison["last_peak_ind_new"] - comparison["last_peak_ind_baseline"]

    comparison["slope_diff"] = comparison["slope_new"] - comparison["slope_baseline"]
    comparison["slope_pct_diff"] = 100 * comparison["slope_diff"] / comparison["slope_baseline"].abs()

    # T-values if present
    if "t_target_pattern_new" in comparison.columns:
        for col in ["t_target_pattern", "t_orth", "t_slope", "d_slope"]:
            comparison[f"{col}_diff"] = comparison[f"{col}_new"] - comparison[f"{col}_baseline"]
            comparison[f"{col}_pct_diff"] = 100 * comparison[f"{col}_diff"] / comparison[f"{col}_baseline"].abs()

    # Summary statistics
    logger.info("\n" + "=" * 70)
    logger.info("COMPARISON WITH BASELINE")
    logger.info("=" * 70)

    logger.info("\nCount_Cells_avg:")
    logger.info(f"  Mean absolute diff: {comparison['Count_Cells_diff'].abs().mean():.4f}")
    logger.info(f"  Mean % diff: {comparison['Count_Cells_pct_diff'].abs().mean():.2f}%")
    logger.info(f"  Max % diff: {comparison['Count_Cells_pct_diff'].abs().max():.2f}%")
    logger.info(f"  Within 1%: {(comparison['Count_Cells_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    logger.info("\nlast_peak_ind:")
    logger.info(f"  Exact matches: {(comparison['last_peak_ind_diff'] == 0).sum()}/{len(comparison)}")
    logger.info(f"  Mean absolute diff: {comparison['last_peak_ind_diff'].abs().mean():.4f}")
    logger.info(f"  Max absolute diff: {comparison['last_peak_ind_diff'].abs().max():.0f}")

    logger.info("\nslope:")
    logger.info(f"  Mean absolute diff: {comparison['slope_diff'].abs().mean():.6f}")
    logger.info(f"  Mean % diff: {comparison['slope_pct_diff'].abs().mean():.2f}%")
    logger.info(f"  Max % diff: {comparison['slope_pct_diff'].abs().max():.2f}%")
    logger.info(f"  Within 10%: {(comparison['slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
    logger.info(f"  Within 1%: {(comparison['slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    # T-value statistics if present
    if "t_target_pattern_new" in comparison.columns:
        logger.info("\nt_target_pattern:")
        logger.info(f"  Mean absolute diff: {comparison['t_target_pattern_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_target_pattern_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_target_pattern_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_target_pattern_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info("\nt_orth:")
        logger.info(f"  Mean absolute diff: {comparison['t_orth_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_orth_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_orth_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_orth_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info("\nt_slope:")
        logger.info(f"  Mean absolute diff: {comparison['t_slope_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['t_slope_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['t_slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['t_slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

        logger.info("\nd_slope:")
        logger.info(f"  Mean absolute diff: {comparison['d_slope_diff'].abs().mean():.6f}")
        logger.info(f"  Mean % diff: {comparison['d_slope_pct_diff'].abs().mean():.2f}%")
        logger.info(f"  Within 10%: {(comparison['d_slope_pct_diff'].abs() < 10).sum()}/{len(comparison)}")
        logger.info(f"  Within 1%: {(comparison['d_slope_pct_diff'].abs() < 1).sum()}/{len(comparison)}")

    logger.info("=" * 70)

    return comparison


def plot_baseline_comparison(dataset: str):
    """Create 2x2 scatter plots (t_target_pattern, slope, t_orth, t_slope).

    Loads comparison CSV and generates PNG with correlation stats.
    Prints greppable correlation summary to logger.
    Requires: comparison CSV from compare_with_baseline_csv()
    """
    logger.info(f"Creating baseline comparison plots for {dataset}")

    # Load comparison CSV
    comparison_path = PROCESSED_DATA_DIR / "virtual_screen_module" / f"{dataset}_baseline_comparison.csv"

    if not comparison_path.exists():
        logger.error(f"Comparison file not found: {comparison_path}")
        logger.error("Run 'haghighi-mito compare-baseline --dataset <dataset>' first")
        return None

    comparison = pd.read_csv(comparison_path)
    n_perts = len(comparison)
    logger.info(f"Loaded {n_perts} perturbations")

    # Check if statistical tests exist
    has_stats = "t_target_pattern_new" in comparison.columns
    if not has_stats:
        logger.error("Statistical tests not found. Run with calculate_stats=True")
        return None

    # Calculate correlations and print greppable output
    metrics = []
    if has_stats:
        for metric in ["t_target_pattern", "t_orth", "t_slope", "slope"]:
            col_new = f"{metric}_new"
            col_baseline = f"{metric}_baseline"

            if col_new in comparison.columns and col_baseline in comparison.columns:
                valid_mask = np.isfinite(comparison[col_new]) & np.isfinite(comparison[col_baseline])
                n_valid = valid_mask.sum()

                if n_valid > 2:
                    corr, _ = pearsonr(comparison.loc[valid_mask, col_baseline], comparison.loc[valid_mask, col_new])

                    pct_diff_col = f"{metric}_pct_diff"
                    if pct_diff_col in comparison.columns:
                        within_10 = (comparison[pct_diff_col].abs() < 10).sum()
                        pct_within_10 = 100 * within_10 / n_perts
                        # GREPPABLE FORMAT - do not change!
                        logger.info(f"{metric}: r={corr:.3f}, within_10%={within_10}/{n_perts} ({pct_within_10:.1f}%)")
                        metrics.append((metric, corr, within_10, pct_within_10))

    # Create scatter plots
    output_dir = PROCESSED_DATA_DIR / "figures" / "diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    def plot_metric(ax, metric_name):
        col_new = f"{metric_name}_new"
        col_baseline = f"{metric_name}_baseline"

        if col_new not in comparison.columns or col_baseline not in comparison.columns:
            ax.text(0.5, 0.5, f"{metric_name} not available", ha="center", va="center", transform=ax.transAxes)
            return

        valid_mask = np.isfinite(comparison[col_new]) & np.isfinite(comparison[col_baseline])
        if valid_mask.sum() < 2:
            ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)
            return

        # Calculate correlation
        corr, _ = pearsonr(comparison.loc[valid_mask, col_baseline], comparison.loc[valid_mask, col_new])

        # Scatter plot
        ax.scatter(comparison.loc[valid_mask, col_baseline], comparison.loc[valid_mask, col_new], alpha=0.5, s=20)

        # Identity line (y=x)
        min_val = min(comparison.loc[valid_mask, col_baseline].min(), comparison.loc[valid_mask, col_new].min())
        max_val = max(comparison.loc[valid_mask, col_baseline].max(), comparison.loc[valid_mask, col_new].max())
        ax.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.5, label="y=x")

        ax.set_xlabel(f"Baseline {metric_name}", fontsize=11)
        ax.set_ylabel(f"Regenerated {metric_name}", fontsize=11)
        ax.set_title(f"{metric_name} (r={corr:.3f})", fontsize=12, fontweight="bold")
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Create all 4 plots
    plot_metric(axes[0, 0], "t_target_pattern")
    plot_metric(axes[0, 1], "slope")
    plot_metric(axes[1, 0], "t_orth")
    plot_metric(axes[1, 1], "t_slope")

    fig.suptitle(f"{dataset.upper()}: Baseline vs Regenerated Metrics (n={n_perts})", fontsize=14, fontweight="bold")
    plt.tight_layout()

    plot_path = output_dir / f"{dataset}_comparison_metrics.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved plots to {plot_path}")
    plt.close()

    return comparison


def compare_with_baseline_csv(dataset: str):
    """Load module + baseline CSVs, compute diffs, save comparison CSV + plots.

    Fast (<1 sec) - generates both CSV and 2x2 scatter plot PNG automatically.
    Prints summary stats, top 5 mismatches, and correlation metrics.
    Raises FileNotFoundError if module CSV missing (run virtual-screen first).
    """
    if dataset not in DATASET_INFO:
        raise ValueError(f"Unknown dataset: {dataset}. Must be one of {list(DATASET_INFO.keys())}")

    logger.info(f"Comparing {dataset} results with baseline...")

    # Load module-generated results
    module_dir = PROCESSED_DATA_DIR / "virtual_screen_module"
    results_path = module_dir / f"{dataset}_results_pattern_aug_070624.csv"

    if not results_path.exists():
        raise FileNotFoundError(f"Module results not found: {results_path}\nRun 'haghighi-mito virtual-screen --dataset {dataset}' first")

    results = pd.read_csv(results_path)
    logger.info(f"Loaded {len(results)} perturbations from {results_path}")

    # Compare with baseline
    comparison = compare_with_baseline(results, dataset)

    # Save comparison
    comparison_path = module_dir / f"{dataset}_baseline_comparison.csv"
    comparison.to_csv(comparison_path, index=False)
    logger.info(f"\nSaved comparison to {comparison_path}")

    # Show some examples of large differences
    has_stats = "t_target_pattern_pct_diff" in comparison.columns
    if has_stats:
        logger.info("\nTop 5 largest t_target_pattern % differences:")
        pert_col = DATASET_INFO[dataset]["pert_col"]
        top_diffs = comparison.nlargest(5, "t_target_pattern_pct_diff")[[pert_col, "t_target_pattern_new", "t_target_pattern_baseline", "t_target_pattern_pct_diff"]]
        print(top_diffs.to_string(index=False))
    else:
        logger.info("\nTop 5 largest slope % differences:")
        pert_col = DATASET_INFO[dataset]["pert_col"]
        top_diffs = comparison.nlargest(5, "slope_pct_diff")[[pert_col, "slope_new", "slope_baseline", "slope_pct_diff"]]
        print(top_diffs.to_string(index=False))

    # Generate diagnostic plots (auto-run)
    logger.info("\nGenerating diagnostic plots...")
    plot_baseline_comparison(dataset)
