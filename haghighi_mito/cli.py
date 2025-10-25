"""Command-line interface for mitochondrial morphology screen analysis."""

from pathlib import Path

import typer
from typing_extensions import Annotated

app = typer.Typer(
    name="haghighi-mito",
    help="Mitochondrial morphology virtual screen analysis pipeline",
)


@app.command(name="process-csvs")
def process_csvs(
    output_dir: Annotated[Path, typer.Option(help="Output directory for Excel files")],
    input_dir: Annotated[Path | None, typer.Option(help="Input directory with CSV files (default: from config)")] = None,
    results_suffix: Annotated[str, typer.Option(help="Suffix for input CSV files")] = "",
    fdr: Annotated[float, typer.Option(help="False discovery rate for BH correction")] = 0.05,
    cell_count_quantile: Annotated[float | None, typer.Option(help="Filter bottom quantile by cell count (e.g., 0.1)")] = None,
):
    """Convert virtual screen CSV results to filtered Excel tables.

    Reads CSV files from virtual screen analysis, applies statistical filtering
    (BH-corrected p-values, orthogonal feature filtering), and saves results
    as multi-sheet Excel files.
    """
    # Lazy imports for faster startup
    from haghighi_mito.config import MITO_VIRTUAL_SCREEN_DIR
    from haghighi_mito.data import convert_virtual_screen_csvs_to_excel

    # Use config default if not provided
    if input_dir is None:
        input_dir = MITO_VIRTUAL_SCREEN_DIR

    # Call function (already has logging via loguru)
    convert_virtual_screen_csvs_to_excel(
        input_dir=input_dir,
        output_dir=output_dir,
        results_suffix=results_suffix,
        fdr=fdr,
        cell_count_quantile=cell_count_quantile,
    )


@app.command(name="create-database")
def create_database_cmd(
    output_path: Annotated[Path, typer.Option(help="Path to output DuckDB database file")],
    tables_dir: Annotated[Path | None, typer.Option(help="Directory containing Excel files (default: from config)")] = None,
    overwrite: Annotated[bool, typer.Option(help="Recreate database even if it exists")] = False,
):
    """Create DuckDB database from Excel screen results.

    Combines all screen result Excel files into a single DuckDB database.
    Only relevant columns are kept for each dataset to reduce database size.
    """
    # Lazy imports for faster startup
    from haghighi_mito.data import create_screen_database

    create_screen_database(
        output_path=output_path,
        tables_dir=tables_dir,
        overwrite=overwrite,
    )


def main():
    """Entry point for CLI."""
    app()


if __name__ == "__main__":
    main()
