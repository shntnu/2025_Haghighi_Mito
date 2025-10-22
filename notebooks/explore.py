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
    # Try local DuckDB file first, fall back to Figshare CSV
    local_db_path = Path(__file__).parent.parent / "data" / "processed" / "screen_results.duckdb"

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

    if local_db_path.exists():
        # Load from local DuckDB
        con = duckdb.connect(str(local_db_path), read_only=True)
        query = f"SELECT {', '.join(common_cols)} FROM screens"
        print(f"Loaded from local DuckDB: {local_db_path}")
    else:
        # Load from Figshare DuckDB using httpfs extension
        figshare_url = "https://figshare.com/ndownloader/files/58970131"
        con = duckdb.connect(':memory:', read_only=False)
        con.execute('INSTALL httpfs')
        con.execute('LOAD httpfs')
        con.execute(f"ATTACH '{figshare_url}' AS remote_db (READ_ONLY)")
        query = f"SELECT {', '.join(common_cols)} FROM remote_db.screens"
        print(f"Loaded from Figshare: {figshare_url}")

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


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Hit calling strategy

    Based on the manuscript methods (sections: "Analysis for virtual screens" and Figure 6):

    ### Primary criteria for calling hits:

    1. **Target phenotype significance** (`p_slope_std`):
       - Must be statistically significant after Benjamini-Hochberg correction
       - Dataset-specific thresholds (critical values):
         - LINCS: 0.00761
         - CDRP: 0.00924
         - JUMP-ORF: 0.00198
         - JUMP-CRISPR: 0.007
         - JUMP-Compound: 0.00365
         - TA-ORF: 0.00409

    2. **Orthogonal features filter** (`p_orth_std`):
       - Must remain above the BH-corrected dataset-specific threshold
       - Ensures non-mitochondrial features stay statistically indistinguishable from controls
       - Avoids compounds that dramatically change overall cell morphology
       - Uses Hotelling's T² test on orthogonal feature set
       - Dataset-specific thresholds (critical values):
         - LINCS: 0.04829
         - CDRP: 0.04812
         - JUMP-ORF: 0.03752
         - JUMP-CRISPR: 0.04785
         - JUMP-Compound: 0.04746
         - TA-ORF: 0.03287

    3. **Effect size** (`d_slope`):
       - Cohen's d of t-test (control vs. perturbation)
       - Positive d_slope: mitochondria more peripheral (depression-like)
       - Negative d_slope: mitochondria more central (psychosis-like)
       - Ranked by absolute d_slope value

    ### Optional filters:
    - Cell count: Avoid perturbations that kill cells
    - Typically filter out bottom 10% by cell count

    ### Definition of orthogonal features:
    - Weak correlation with MITO-SLOPE (|r| < 0.3)
    - Not distinctive across patient groups
    - Low inter-feature redundancy (|r| < 0.8)

    ### Hit categories:
    - **Positive hits**: d_slope > 0, potential depression phenotype modulators
    - **Negative hits**: d_slope < 0, potential psychosis phenotype modulators
    """
    )
    return


@app.cell
def _(filtered_df, pl):
    # Dataset-specific BH-corrected thresholds for target and orthogonal tests
    # Keys must match Metadata_dataset values in screens table
    target_bh_corrected_critical_dict = {
        'LINCS': 0.00761,
        'CDRP': 0.00924,
        'JUMP_ORF': 0.00198,
        'JUMP_CRISPR': 0.007,
        'JUMP_Compound': 0.00365,
        'TA_ORF': 0.00409
    }
    orth_bh_corrected_critical_dict = {
        'LINCS': 0.04829,
        'CDRP': 0.04812,
        'JUMP_ORF': 0.03752,
        'JUMP_CRISPR': 0.04785,
        'JUMP_Compound': 0.04746,
        'TA_ORF': 0.03287
    }

    # Apply hit calling criteria with dataset-specific thresholds
    hits = filtered_df.filter(
        # Filter 1: Orthogonal features must remain normal (p_orth_std > threshold)
        # Filter 2: Target phenotype must be significant (p_slope_std < threshold)
        (
            (pl.col('Metadata_dataset') == 'LINCS') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['LINCS']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['LINCS'])
        ) |
        (
            (pl.col('Metadata_dataset') == 'CDRP') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['CDRP']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['CDRP'])
        ) |
        (
            (pl.col('Metadata_dataset') == 'JUMP_ORF') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['JUMP_ORF']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['JUMP_ORF'])
        ) |
        (
            (pl.col('Metadata_dataset') == 'JUMP_CRISPR') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['JUMP_CRISPR']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['JUMP_CRISPR'])
        ) |
        (
            (pl.col('Metadata_dataset') == 'JUMP_Compound') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['JUMP_Compound']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['JUMP_Compound'])
        ) |
        (
            (pl.col('Metadata_dataset') == 'TA_ORF') &
            (pl.col('p_orth_std') > orth_bh_corrected_critical_dict['TA_ORF']) &
            (pl.col('p_slope_std') < target_bh_corrected_critical_dict['TA_ORF'])
        )
    ).sort('d_slope')  # Sort by effect size

    # Separate into positive and negative hits
    negative_hits = hits.filter(pl.col('d_slope') < 0)  # Psychosis-like (mitochondria more central)
    positive_hits = hits.filter(pl.col('d_slope') > 0)  # Depression-like (mitochondria more peripheral)
    return hits, negative_hits, positive_hits


@app.cell
def _(hits, mo, pl):
    mo.md(
        f"""
    ## Hits Summary

    Total hits meeting **both** criteria: **{len(hits):,}**
    - Negative hits (psychosis-like): {len(hits.filter(pl.col('d_slope') < 0)):,}
    - Positive hits (depression-like): {len(hits.filter(pl.col('d_slope') > 0)):,}
    """
    )
    return


@app.cell
def _(mo):
    mo.md("""#### Hits by Dataset""")
    return


@app.cell
def _(hits, mo, pl):
    hits_by_dataset = hits.group_by('Metadata_dataset').agg([
        pl.len().alias('total_hits'),
        pl.col('d_slope').filter(pl.col('d_slope') < 0).len().alias('negative_hits'),
        pl.col('d_slope').filter(pl.col('d_slope') > 0).len().alias('positive_hits')
    ]).sort('total_hits', descending=True)
    mo.ui.table(hits_by_dataset, selection=None)
    return (hits_by_dataset,)


@app.cell
def _(mo):
    mo.md(
        """
    ### Negative Hits (Psychosis-like: d_slope < 0)

    Mitochondria more centralized, similar to psychosis patient phenotype
    """
    )
    return


@app.cell
def _(mo, negative_hits):
    mo.ui.table(negative_hits, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        'p_orth_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


@app.cell
def _(mo):
    mo.md(
        """
    ### Positive Hits (Depression-like: d_slope > 0)

    Mitochondria more peripheral, similar to depression patient phenotype
    """
    )
    return


@app.cell
def _(mo, positive_hits):
    mo.ui.table(positive_hits, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        'p_orth_std': '{:.6f}'.format,
        't_orth': '{:.2f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


@app.cell
def _(Path, hits, hits_by_dataset, negative_hits, positive_hits):
    # Save results to CSV only if directory exists
    output_dir = Path(__file__).parent.parent / "data" / "processed"

    if output_dir.exists():
        hits_by_dataset.write_csv(output_dir / "hits_by_dataset.csv")
        hits.write_csv(output_dir / "all_hits.csv")
        negative_hits.write_csv(output_dir / "negative_hits.csv")
        positive_hits.write_csv(output_dir / "positive_hits.csv")
        print(f"Saved hits to {output_dir}")
    else:
        print(f"Output directory does not exist: {output_dir}")
    return


@app.cell
def _(mo):
    mo.md("""## OTC/GRAS Compounds Analysis""")
    return


@app.cell
def _(mo):
    mo.md(
        """
    Searching for over-the-counter (OTC) medications and Generally Recognized As Safe (GRAS)
    compounds/supplements in the LINCS dataset to identify readily available compounds with
    mitochondrial localization effects.
    """
    )
    return


@app.cell
def _(pl, screens):
    # LINCS dataset-specific BH-corrected thresholds (FDR=0.05)
    LINCS_TARGET_THRESHOLD = 0.00761  # Target phenotype must be significant
    LINCS_ORTH_THRESHOLD = 0.04829    # Orthogonal features must remain normal

    # List of OTC/GRAS compounds to search for
    otc_gras_compounds = [
        'caffeine', 'melatonin', 'aspirin', 'ibuprofen', 'naproxen',
        'metformin', 'niacin', 'resveratrol', 'curcumin', 'quercetin',
        'ranitidine', 'omeprazole', 'cimetidine', 'famotidine',
        'paracetamol', 'dextromethorphan', 'loperamide', 'loratadine',
        'levocetirizine', 'phenylephrine', 'cyanocobalamin', 'folic-acid',
        'ginkgolide-b'
    ]

    # Filter for LINCS dataset and OTC/GRAS compounds
    otc_results = screens.filter(
        (pl.col('Metadata_dataset') == 'LINCS') &
        (pl.col('Metadata_pert_id').is_in(otc_gras_compounds))
    ).with_columns([
        # Add filter status column
        pl.when(
            (pl.col('p_slope_std') < LINCS_TARGET_THRESHOLD) &
            (pl.col('p_orth_std') > LINCS_ORTH_THRESHOLD)
        ).then(pl.lit('PASS'))
        .when(pl.col('p_orth_std') > LINCS_ORTH_THRESHOLD)
        .then(pl.lit('Near miss (p_slope)'))
        .when(pl.col('p_slope_std') < LINCS_TARGET_THRESHOLD)
        .then(pl.lit('Near miss (p_orth)'))
        .otherwise(pl.lit('Fail both'))
        .alias('filter_status')
    ]).select([
        'Metadata_pert_id',
        'd_slope',
        'p_slope_std',
        'p_orth_std',
        'Count_Cells_avg',
        'filter_status'
    ]).sort(['filter_status', 'd_slope'], descending=[False, True])
    return (otc_results,)


@app.cell
def _(mo, otc_results, pl):
    passing = otc_results.filter(pl.col('filter_status') == 'PASS')
    near_miss = otc_results.filter(pl.col('filter_status').str.contains('Near miss'))

    mo.md(
        f"""
    ### Results Summary

    - **Passing filters**: {len(passing)}
    - **Near misses**: {len(near_miss)}
    - **Total tested**: {len(otc_results)}
    """
    )
    return


@app.cell
def _(mo, otc_results):
    mo.ui.table(otc_results, selection=None, format_mapping={
        'd_slope': '{:.4f}'.format,
        'p_slope_std': '{:.6f}'.format,
        'p_orth_std': '{:.6f}'.format,
        'Count_Cells_avg': '{:.2f}'.format
    })
    return


if __name__ == "__main__":
    app.run()
