import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import duckdb
    import polars as pl
    from pathlib import Path
    return Path, duckdb, mo, pl


@app.cell
def _(mo):
    mo.md("""# Screen Results Explorer""")
    return


@app.cell
def _(Path, duckdb, pl):
    # Use absolute path relative to this notebook
    db_path = Path(__file__).parent.parent / 'data' / 'processed' / 'screen_results.duckdb'
    con = duckdb.connect(str(db_path), read_only=True)

    # Load only common columns for faster, cleaner analysis
    common_cols = [
        'Metadata_dataset',
        'Metadata_pert_type',
        'Metadata_pert_id',
        'd_slope',
        'p_slope_std',
        'p_orth_std',
        't_orth',
        'Count_Cells_avg'
    ]

    query = f"SELECT {', '.join(common_cols)} FROM screens"
    screens = pl.from_arrow(con.execute(query).fetch_arrow_table())
    return (screens,)


@app.cell
def _(screens):
    summary_stats = {
        'total_records': len(screens),
        'n_datasets': screens['Metadata_dataset'].n_unique()
    }
    return (summary_stats,)


@app.cell
def _(mo, summary_stats):
    mo.md(
        f"""
    Total records: {summary_stats['total_records']:,}

    Datasets: {summary_stats['n_datasets']}
    """
    )
    return


@app.cell
def _(mo, screens):
    # Stratified random sample: 2 from each dataset
    sample = screens.group_by('Metadata_dataset').map_groups(
        lambda group: group.sample(n=min(2, len(group)), seed=42)
    ).sort('Metadata_dataset')
    mo.ui.table(sample, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        'p_orth_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


@app.cell
def _(mo):
    mo.md("""## Dataset Breakdown""")
    return


@app.cell
def _(mo, screens):
    dataset_counts = screens['Metadata_dataset'].value_counts().sort('count', descending=True)
    mo.ui.table(dataset_counts, selection=None, format_mapping={'count': '{:,}'.format})
    return


@app.cell
def _(mo):
    mo.md("""## Key Metrics""")
    return


@app.cell
def _(screens):
    screens.select(['d_slope', 'p_slope_std', 'Count_Cells_avg']).describe()
    return


@app.cell
def _(mo):
    mo.md("""## Filter by Dataset""")
    return


@app.cell
def _(mo, screens):
    dataset_list = sorted(screens['Metadata_dataset'].drop_nulls().unique().to_list())

    dataset_dropdown = mo.ui.dropdown(
        options=['All'] + dataset_list,
        value='All',
        label='Select Dataset:'
    )
    dataset_dropdown
    return (dataset_dropdown,)


@app.cell
def _(dataset_dropdown, pl, screens):
    if dataset_dropdown.value == 'All':
        filtered_df = screens
    else:
        filtered_df = screens.filter(pl.col('Metadata_dataset') == dataset_dropdown.value)
    return (filtered_df,)


@app.cell
def _(filtered_df, mo, pl):
    null_count = filtered_df.filter(pl.col('d_slope').is_null()).height
    mo.md(
        f"""
    ### Filtered Results
    Showing {len(filtered_df):,} records

    {"⚠️ **Warning:** " + str(null_count) + " records with null d_slope values" if null_count > 0 else ""}
    """
    )
    return (null_count,)


@app.cell
def _(mo, null_count):
    if null_count > 0:
        mo.md("""### Data Quality Issues - Null d_slope Values""")
    return


@app.cell
def _(filtered_df, mo, null_count, pl):
    if null_count > 0:
        null_records = filtered_df.filter(pl.col('d_slope').is_null())
        mo.ui.table(null_records, selection=None)
    return


@app.cell
def _(mo):
    mo.md("""## Top Hits""")
    return


@app.cell
def _(mo):
    mo.md("""**Most negative d_slope** (mitochondria more centralized, psychosis-like)""")
    return


@app.cell
def _(filtered_df, mo, pl):
    # Filter out nulls before sorting
    top_negative = filtered_df.filter(
        pl.col('d_slope').is_not_null()
    ).sort('d_slope').head(10)

    mo.ui.table(top_negative, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


@app.cell
def _(mo):
    mo.md("""**Most positive d_slope** (mitochondria more peripheral, depression-like)""")
    return


@app.cell
def _(filtered_df, mo, pl):
    # Filter out nulls before sorting
    top_positive = filtered_df.filter(
        pl.col('d_slope').is_not_null()
    ).sort('d_slope', descending=True).head(10)

    mo.ui.table(top_positive, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


@app.cell
def _(mo):
    mo.md("""## Significant Hits (p < 0.05)""")
    return


@app.cell
def _(filtered_df, mo, pl):
    sig_hits = filtered_df.filter(pl.col('p_slope_std') < 0.05)
    mo.ui.table(sig_hits.head(100), selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return (sig_hits,)


@app.cell
def _(mo, sig_hits):
    mo.md(f"""Found **{len(sig_hits):,}** significant hits (p_slope_std < 0.05)""")
    return


if __name__ == "__main__":
    app.run()
