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
            virtual_screen_module/{dataset}_comparison_metrics.png (plots)
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
    Also preserves provenance metadata (n_sites, n_plates, etc.) from regenerated data.
    """
    baseline_path = EXTERNAL_DATA_DIR / "mito_project/workspace/results/virtual_screen_baseline" / f"{dataset}_results_pattern_aug_070624.csv"

    logger.info(f"Loading baseline from {baseline_path}")
    baseline = pd.read_csv(baseline_path)

    pert_col = DATASET_INFO[dataset]["pert_col"]
    baseline_cols = [pert_col, "Count_Cells_avg", "last_peak_ind", "slope"]

    # Add t-values if they exist
    if "t_target_pattern" in results.columns:
        baseline_cols.extend(["t_target_pattern", "t_orth", "t_slope", "d_slope"])

    # Note: Provenance columns (n_sites, n_plates, etc.) are NOT in baseline
    # They will be kept from the 'results' dataframe (regenerated data) with no suffix
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
    output_dir = PROCESSED_DATA_DIR / "virtual_screen_module"
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


def _compute_summary(comparison: pd.DataFrame) -> pd.DataFrame:
    """Compute summary statistics for each metric.

    Also computes provenance summary if available (n_sites, n_plates, etc.).
    """
    metrics = ["Count_Cells_avg", "slope", "last_peak_ind"]
    if "t_target_pattern_new" in comparison.columns:
        metrics.extend(["t_target_pattern", "t_orth", "t_slope", "d_slope"])

    rows = []
    for metric in metrics:
        col_new = f"{metric}_new"
        col_baseline = f"{metric}_baseline"

        if col_new not in comparison.columns or col_baseline not in comparison.columns:
            continue

        valid_mask = np.isfinite(comparison[col_new]) & np.isfinite(comparison[col_baseline])
        n_valid = valid_mask.sum()

        if n_valid < 2:
            continue

        # Correlations
        pearson_r, _ = pearsonr(comparison.loc[valid_mask, col_baseline], comparison.loc[valid_mask, col_new])
        spearman_r, _ = spearmanr(comparison.loc[valid_mask, col_baseline], comparison.loc[valid_mask, col_new])

        # Differences (handle special case for Count_Cells_avg)
        if metric == "Count_Cells_avg":
            diff_col = "Count_Cells_diff"
            pct_diff_col = "Count_Cells_pct_diff"
        else:
            diff_col = f"{metric}_diff"
            pct_diff_col = f"{metric}_pct_diff"

        # Use absolute differences from the _diff column
        if diff_col in comparison.columns:
            mean_abs_diff = comparison[diff_col].abs().mean()
            median_abs_diff = comparison[diff_col].abs().median()
            max_abs_diff = comparison[diff_col].abs().max()
        else:
            mean_abs_diff = median_abs_diff = max_abs_diff = np.nan

        # Within thresholds (for pct_diff where applicable)
        if pct_diff_col in comparison.columns:
            within_1pct = (comparison[pct_diff_col].abs() < 1).sum()
            within_10pct = (comparison[pct_diff_col].abs() < 10).sum()
            pct_within_10pct = 100 * within_10pct / len(comparison)
        else:
            within_1pct = within_10pct = pct_within_10pct = np.nan

        rows.append({
            "metric": metric,
            "n_samples": len(comparison),
            "pearson_r": pearson_r,
            "spearman_r": spearman_r,
            "mean_abs_diff": mean_abs_diff,
            "median_abs_diff": median_abs_diff,
            "max_abs_diff": max_abs_diff,
            "within_1pct": within_1pct,
            "within_10pct": within_10pct,
            "pct_within_10pct": pct_within_10pct,
        })

    return pd.DataFrame(rows)


def _compute_provenance_summary(comparison: pd.DataFrame) -> dict:
    """Compute summary statistics for provenance metadata.

    Returns dict with summary stats for n_sites, n_plates, n_wells, etc.
    Returns empty dict if provenance columns not available.
    """
    provenance_cols = ["n_sites", "n_plates", "n_wells", "Count_Cells_std", "slope_std"]

    # Check if any provenance columns exist
    available_cols = [col for col in provenance_cols if col in comparison.columns]

    if not available_cols:
        return {}

    summary = {}
    for col in available_cols:
        if col in comparison.columns:
            summary[f"{col}_mean"] = comparison[col].mean()
            summary[f"{col}_median"] = comparison[col].median()
            summary[f"{col}_min"] = comparison[col].min()
            summary[f"{col}_max"] = comparison[col].max()

    # Compute low sample size warnings
    if "n_sites" in comparison.columns:
        n_low_sites = (comparison["n_sites"] < 10).sum()
        summary["n_low_sites_count"] = n_low_sites
        summary["n_low_sites_pct"] = 100 * n_low_sites / len(comparison)

    if "n_plates" in comparison.columns:
        n_single_plate = (comparison["n_plates"] == 1).sum()
        summary["n_single_plate_count"] = n_single_plate
        summary["n_single_plate_pct"] = 100 * n_single_plate / len(comparison)

    return summary


def _format_summary(summary: pd.DataFrame, dataset: str, provenance: dict = None) -> str:
    """Format summary statistics as a readable table.

    Args:
        summary: Metric comparison statistics DataFrame
        dataset: Dataset name
        provenance: Optional dict with provenance summary stats
    """
    lines = []
    lines.append(f"{dataset.upper()}: Baseline Comparison Summary (n={summary['n_samples'].iloc[0] if len(summary) > 0 else 0})")
    lines.append("=" * 80)
    lines.append(f"{'Metric':<20} {'Pearson r':>10} {'Within 10%':>12} {'Mean Abs Diff':>15}")
    lines.append("-" * 80)

    for _, row in summary.iterrows():
        metric = row["metric"]
        r = row["pearson_r"]
        within_10 = row["pct_within_10pct"]
        mean_diff = row["mean_abs_diff"]

        # Add warning indicators
        warn = ""
        if pd.notna(within_10):
            if within_10 < 20:
                warn = "  ⚠⚠"
            elif within_10 < 50:
                warn = "  ⚠"

        within_str = f"{within_10:.1f}%" if pd.notna(within_10) else "N/A"
        lines.append(f"{metric:<20} {r:>10.3f} {within_str:>12} {mean_diff:>15.3f}{warn}")

    lines.append("=" * 80)
    lines.append("⚠ = <50% within 10%, ⚠⚠ = <20% within 10%")

    # Add provenance section if available
    if provenance and len(provenance) > 0:
        lines.append("")
        lines.append("PROVENANCE METADATA (regenerated data only)")
        lines.append("-" * 80)

        if "n_sites_median" in provenance:
            lines.append(f"Sites per perturbation:  median={provenance['n_sites_median']:.0f}, "
                        f"range={provenance['n_sites_min']:.0f}-{provenance['n_sites_max']:.0f}")
            if "n_low_sites_count" in provenance:
                lines.append(f"  ⚠ {provenance['n_low_sites_count']:.0f} perturbations ({provenance['n_low_sites_pct']:.1f}%) have <10 sites")

        if "n_plates_median" in provenance:
            lines.append(f"Plates per perturbation: median={provenance['n_plates_median']:.0f}, "
                        f"range={provenance['n_plates_min']:.0f}-{provenance['n_plates_max']:.0f}")
            if "n_single_plate_count" in provenance:
                lines.append(f"  ⚠ {provenance['n_single_plate_count']:.0f} perturbations ({provenance['n_single_plate_pct']:.1f}%) have only 1 plate")

        if "slope_std_median" in provenance:
            lines.append(f"Slope variability:       median={provenance['slope_std_median']:.3f}, "
                        f"range={provenance['slope_std_min']:.3f}-{provenance['slope_std_max']:.3f}")

    return "\n".join(lines)


def export_examples(comparison: pd.DataFrame, dataset: str, metric: str = "slope", n_per_tail: int = 10):
    """Export best and worst matches (extremes) with provenance for comparison.

    Exports both tails of the distribution:
    - Top n_per_tail worst matches (largest absolute % difference)
    - Top n_per_tail best matches (smallest absolute % difference)

    This allows comparing what distinguishes good vs poor reproducibility.

    Args:
        comparison: Comparison DataFrame with *_new, *_baseline, *_diff columns
        dataset: Dataset name
        metric: Which metric to analyze (default: slope, the most important)
        n_per_tail: Number of examples per tail to export (default: 10)

    Saves: {dataset}_{metric}_examples.csv with best/worst matches + provenance
    """
    module_dir = PROCESSED_DATA_DIR / "virtual_screen_module"

    # Get columns for this metric
    new_col = f"{metric}_new"
    baseline_col = f"{metric}_baseline"
    diff_col = f"{metric}_diff"
    pct_diff_col = f"{metric}_pct_diff"

    # Check columns exist
    if new_col not in comparison.columns or baseline_col not in comparison.columns:
        logger.warning(f"Cannot export examples: {new_col} or {baseline_col} not in comparison")
        return

    # Calculate absolute percentage difference and add as temporary column
    if pct_diff_col in comparison.columns:
        comparison["_abs_pct_diff"] = comparison[pct_diff_col].abs()
        sort_col = "_abs_pct_diff"
    else:
        comparison["_abs_diff"] = comparison[diff_col].abs()
        sort_col = "_abs_diff"

    # Get worst matches (largest absolute % difference)
    worst = comparison.nlargest(n_per_tail, sort_col).copy()
    worst["match_quality"] = "worst"

    # Get best matches (smallest absolute % difference)
    best = comparison.nsmallest(n_per_tail, sort_col).copy()
    best["match_quality"] = "best"

    # Remove temporary column
    comparison.drop(columns=[sort_col], inplace=True)

    # Combine and sort by match quality, then by difference
    examples = pd.concat([worst, best], ignore_index=True)

    # Select relevant columns for output
    meta_cols = DATASET_INFO[dataset]["meta_cols"][:2]  # Just first 2 metadata cols (gene_name, pert_name)

    output_cols = ["match_quality"] + meta_cols.copy()
    output_cols.extend([baseline_col, new_col, diff_col])

    if pct_diff_col in comparison.columns:
        output_cols.append(pct_diff_col)

    # Add provenance columns if available
    provenance_cols = ["n_sites", "n_plates", "n_wells", "Count_Cells_std", "slope_std", "median_plate_id"]
    for col in provenance_cols:
        if col in examples.columns:
            output_cols.append(col)

    # Filter to existing columns and export
    output_cols = [col for col in output_cols if col in examples.columns]
    examples_export = examples[output_cols].copy()

    # Round numeric columns to 6 decimal places for readability
    numeric_cols = examples_export.select_dtypes(include=[np.number]).columns
    examples_export[numeric_cols] = examples_export[numeric_cols].round(6)

    output_path = module_dir / f"{dataset}_{metric}_examples.csv"
    examples_export.to_csv(output_path, index=False)
    logger.info(f"Exported {n_per_tail} best + {n_per_tail} worst {metric} matches to {output_path}")


def compare_with_baseline_csv(dataset: str):
    """Load module + baseline CSVs, compute diffs, save comparison CSV + plots.

    Fast (<1 sec) - generates both CSV and 2x2 scatter plot PNG automatically.
    Saves diagnostic summary (CSV + TXT) and correlation plots.
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
    logger.info(f"Saved comparison to {comparison_path}")

    # Generate summary statistics
    summary = _compute_summary(comparison)
    provenance_summary = _compute_provenance_summary(comparison)

    summary_csv_path = module_dir / f"{dataset}_diagnostic_summary.csv"
    summary_txt_path = module_dir / f"{dataset}_diagnostic_summary.txt"

    # Round numeric columns to 6 decimal places for readability
    summary_rounded = summary.round(6)
    summary_rounded.to_csv(summary_csv_path, index=False)

    summary_txt = _format_summary(summary, dataset, provenance_summary)
    summary_txt_path.write_text(summary_txt)

    logger.info(f"Saved diagnostic summary to {summary_csv_path}")
    logger.info(f"Saved formatted summary to {summary_txt_path}")

    # Generate diagnostic plots
    plot_baseline_comparison(dataset)

    # Export examples (best + worst matches) for slope analysis
    export_examples(comparison, dataset, metric="slope", n_per_tail=10)

    # Print summary to console
    print("\n" + summary_txt)
