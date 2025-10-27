import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import altair as alt
    from pathlib import Path
    return Path, alt, mo, pl


@app.cell
def _(mo):
    mo.md("""# Baseline vs Regenerated Comparison""")
    return


@app.cell
def _(Path):
    comparison_dir = Path("/Users/shsingh/Documents/GitHub/misc/2025_Haghighi_Mito/data/processed/virtual_screen_module")
    dataset_options = [f.stem.replace("_baseline_comparison", "") for f in comparison_dir.glob("*_baseline_comparison.csv")]
    return comparison_dir, dataset_options


@app.cell
def _(dataset_options, mo):
    dataset_selector = mo.ui.dropdown(options=dataset_options, value=dataset_options[0] if dataset_options else None, label="Dataset")
    metric_selector = mo.ui.dropdown(options=["slope", "t_target_pattern", "t_orth", "t_slope", "d_slope"], value="slope", label="Metric")
    mo.hstack([dataset_selector, metric_selector])
    return dataset_selector, metric_selector


@app.cell
def _(alt, comparison_dir, dataset_selector, metric_selector, mo, pl):
    if not dataset_selector.value:
        _output = mo.md("Please select a dataset")
    else:
        _df = pl.read_csv(comparison_dir / f"{dataset_selector.value}_baseline_comparison.csv")
        _baseline_col = f"{metric_selector.value}_baseline"
        _new_col = f"{metric_selector.value}_new"

        _df = _df.filter(
            pl.col(_baseline_col).is_not_null() & pl.col(_new_col).is_not_null() &
            pl.col(_baseline_col).is_finite() & pl.col(_new_col).is_finite()
        )

        _min_val = min(_df[_baseline_col].min(), _df[_new_col].min())
        _max_val = max(_df[_baseline_col].max(), _df[_new_col].max())

        _scatter = alt.Chart(_df).mark_circle(opacity=0.5, size=30).encode(
            x=alt.X(_baseline_col, title=f"Baseline {metric_selector.value}", scale=alt.Scale(domain=[_min_val, _max_val])),
            y=alt.Y(_new_col, title=f"Regenerated {metric_selector.value}", scale=alt.Scale(domain=[_min_val, _max_val])),
            tooltip=["Metadata_gene_name:N", alt.Tooltip(_baseline_col, format=".4f"), alt.Tooltip(_new_col, format=".4f"), alt.Tooltip(f"{metric_selector.value}_pct_diff", format=".2f")]
        )

        _identity = alt.Chart(pl.DataFrame({"x": [_min_val, _max_val], "y": [_min_val, _max_val]})).mark_line(color='red', strokeDash=[5,5]).encode(x='x:Q', y='y:Q')

        _chart = (_scatter + _identity).properties(title=f"{dataset_selector.value}: {metric_selector.value}", width=500, height=500).interactive()

        _n = len(_df)
        _corr = _df.select(pl.corr(_baseline_col, _new_col)).item()
        _within_10 = _df.filter(pl.col(f"{metric_selector.value}_pct_diff").abs() <= 10)
        _match_rate = len(_within_10) / _n * 100

        _output = mo.vstack([_chart, mo.md(f"**N = {_n} | Correlation = {_corr:.4f} | Within 10% = {_match_rate:.1f}% ({len(_within_10)}/{_n})**")])

    _output
    return


if __name__ == "__main__":
    app.run()
