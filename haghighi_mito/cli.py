"""Command-line interface for mitochondrial morphology screen analysis."""

from pathlib import Path

import typer
from typing_extensions import Annotated

app = typer.Typer(
    name="haghighi-mito",
    help="Mitochondrial morphology virtual screen analysis pipeline",
)


@app.command(name="process-csv-single")
def process_csv_single(
    dataset: Annotated[str, typer.Option(help="Dataset name (e.g., 'CDRP', 'lincs', 'jump_orf')")],
    csv_path: Annotated[Path, typer.Option(help="Path to input CSV file")],
    output_dir: Annotated[Path, typer.Option(help="Output directory for Excel file")],
    fdr: Annotated[float, typer.Option(help="False discovery rate for BH correction")] = 0.05,
    cell_count_quantile: Annotated[float | None, typer.Option(help="Filter bottom quantile by cell count (e.g., 0.1)")] = None,
    parquet_output_dir: Annotated[Path | None, typer.Option(help="Directory to save Parquet file (optional)")] = None,
):
    """Process a single virtual screen CSV file to Excel and Parquet.

    Reads a single CSV file from virtual screen analysis, applies statistical filtering
    (BH-corrected p-values, orthogonal feature filtering), and saves results
    as a multi-sheet Excel file. Optionally saves unfiltered data as Parquet.
    """
    # Lazy imports for faster startup
    from haghighi_mito.data import process_single_virtual_screen_csv

    # Call function (already has logging via loguru)
    process_single_virtual_screen_csv(
        dataset=dataset,
        csv_path=csv_path,
        output_dir=output_dir,
        fdr=fdr,
        cell_count_quantile=cell_count_quantile,
        parquet_output_dir=parquet_output_dir,
    )


@app.command(name="create-database")
def create_database_cmd(
    output_path: Annotated[Path, typer.Option(help="Path to output DuckDB database file")],
    tables_dir: Annotated[Path | None, typer.Option(help="Directory containing Excel files (default: from config)")] = None,
    overwrite: Annotated[bool, typer.Option(help="Recreate database even if it exists")] = False,
    use_parquet: Annotated[bool, typer.Option(help="Read from Parquet files instead of Excel")] = False,
    parquet_dir: Annotated[Path | None, typer.Option(help="Directory containing Parquet files (required if --use-parquet)")] = None,
):
    """Create DuckDB database from Excel or Parquet screen results.

    Combines all screen result files into a single DuckDB database.
    Only relevant columns are kept for each dataset to reduce database size.
    """
    # Lazy imports for faster startup
    from haghighi_mito.data import create_screen_database

    create_screen_database(
        output_path=output_path,
        tables_dir=tables_dir,
        overwrite=overwrite,
        use_parquet=use_parquet,
        parquet_dir=parquet_dir,
    )


@app.command(name="validate-databases")
def validate_databases_cmd(
    baseline: Annotated[Path, typer.Option(help="Path to baseline DuckDB file")],
    new: Annotated[Path, typer.Option(help="Path to new DuckDB file to validate")],
):
    """Compare two DuckDB databases for identical data.

    Validates that two DuckDB databases contain identical data by comparing
    all rows and columns. Useful for validating refactoring or pipeline changes.
    """
    import sys

    # Lazy imports for faster startup
    from haghighi_mito.data import validate_databases

    identical = validate_databases(baseline_path=baseline, new_path=new)

    if not identical:
        sys.exit(1)


@app.command(name="virtual-screen")
def virtual_screen_cmd(
    dataset: Annotated[str, typer.Option(help="Dataset to analyze (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)")],
    compare_baseline: Annotated[bool, typer.Option(help="Compare results with baseline CSV")] = True,
    calculate_stats: Annotated[bool, typer.Option(help="Calculate statistical tests (t-values)")] = True,
):
    """Run virtual screen analysis from scratch.

    Calculates radial distribution metrics (Count_Cells_avg, last_peak_ind, slope)
    and statistical tests (t_target_pattern, t_orth, t_slope, d_slope)
    from per-site profiles.
    """
    # Lazy imports for faster startup
    from haghighi_mito.virtual_screen import run_virtual_screen

    run_virtual_screen(dataset=dataset, compare_baseline=compare_baseline, calculate_stats=calculate_stats)


@app.command(name="analyze-t-target-pattern")
def analyze_t_target_pattern_cmd(
    dataset: Annotated[str, typer.Option(help="Dataset to analyze (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)")],
):
    """Analyze t_target_pattern distribution and relationship to baseline.

    Explores the distribution of t_target_pattern (Hotelling's T² on full radial
    distribution) and its relationship to the baseline values. Analyzes:
    - Distribution statistics (mean, median, percentiles)
    - Correlation with baseline (Pearson, Spearman)
    - Systematic transformations (linear regression)
    - Relationship between t_target_pattern match quality and slope match quality
    - Divergence patterns (where t_target_pattern matches but slope diverges)

    t_target_pattern is the most direct test for pattern similarity as it bypasses
    peak detection entirely and compares the full 12-bin radial distribution using
    multivariate statistics.

    Requires running 'virtual-screen --compare-baseline' first to generate comparison file.
    """
    # Lazy imports for faster startup
    from haghighi_mito.virtual_screen import analyze_t_target_pattern_distribution

    analyze_t_target_pattern_distribution(dataset=dataset)


@app.command(name="compare-per-plate")
def compare_per_plate_cmd(
    dataset: Annotated[str, typer.Option(help="Dataset to analyze (taorf, CDRP, lincs, jump_orf, jump_crispr, jump_compound)")],
):
    """Compare per-plate statistical results with baseline.

    Instead of comparing median-aggregated results (one value per perturbation),
    this compares ALL plate-level results. For each perturbation that appears on
    multiple plates, checks if ANY of the regenerated plate values match the
    baseline plate values.

    This tests whether:
    - Per-plate calculations are identical (problem is aggregation method)
    - Or per-plate calculations differ (problem is in core statistics)

    Requires running 'virtual-screen --compare-baseline' first.
    """
    # Lazy imports for faster startup
    from haghighi_mito.virtual_screen import compare_per_plate_results

    compare_per_plate_results(dataset=dataset)


def main():
    """Entry point for CLI."""
    app()


if __name__ == "__main__":
    main()
