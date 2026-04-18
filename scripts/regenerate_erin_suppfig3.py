"""Regenerate Erin's upstream Supp Fig 3 against the corrected label CSV.

Mirrors the structure of Erin's upstream
``plot_singlecell_cellsize.ipynb`` (at upstream commit ``ef8f26b``) but with
two substitutions:

* ``ag_df`` is read from the **corrected** aggregated CSV
  (``data/interim/upstream_corrected/aggregated_profiles_fibroblast.csv``),
  which has ``272→Control`` and ``MCL004→SZA`` applied. Erin's original
  notebook does not apply these overrides.
* The per-subject SQL loop (cell 3) is replaced by a single ``read_parquet``
  against our cached ``data/interim/fibroblast_singlecell.parquet`` (built by
  ``haghighi-mito build-singlecell-parquet``). The parquet is an exact
  concatenation of Cells+Cytoplasm+Nuclei joins from the same 185 SQLite
  backends Erin loops, so downstream rows are identical.

All QC constants, feature list, and plotting style match Erin's notebook
verbatim so the regenerated PDFs can be handed to Erin/Anne as drop-in
replacements for the published Supp Fig 3 panels.

Run: ``.pixi/envs/default/bin/python scripts/regenerate_erin_suppfig3.py``

Outputs: 5 PDFs under ``data/interim/upstream_corrected/erin_suppfig3_corrected/``.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

from haghighi_mito.config import FIBROBLAST_SINGLECELL_PARQUET

# --- inputs ------------------------------------------------------------------
CORRECTED_AGG_CSV = Path("data/interim/upstream_corrected/aggregated_profiles_fibroblast.csv")
OUT_DIR = Path("data/interim/upstream_corrected/erin_suppfig3_corrected")

# --- Erin's cell 17 constants (µm/px and style) ------------------------------
MICRONS_PER_PIXEL = 0.161
PATIENT_ORDER = ["Control", "psychosis", "BP", "SZ", "SZA", "MDD"]
PALETTE = ["lightgray", "#a484ac", "firebrick", "lightcoral", "pink", "lightsteelblue"]

# --- Erin's cell 5 QC constants ----------------------------------------------
BORDER_PX = 200
IM_WIDTH = 1388
IM_HEIGHT = 1040
NUCLEI_MAJOR_RATIO_MIN = 2
NUCLEI_AREA_RATIO_MIN = 5
CELLS_ACTIN_MAX = 0.5
NUCLEI_ACTIN_MAX = 0.55


def round_sig(x: float, sig: int = 2) -> float:
    """Mirror Erin's ``round_sig`` used for annotated p-values."""
    import math

    if x == 0 or not math.isfinite(x):
        return x
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sig - 1))


def load_corrected_ag_df() -> pd.DataFrame:
    """Cell 2 analogue: read the corrected aggregated CSV as the label source."""
    ag = pd.read_csv(CORRECTED_AGG_CSV)
    ag["subject"] = ag["subject"].astype(str)
    return ag


def load_raw_single_cells(ag_df: pd.DataFrame) -> pd.DataFrame:
    """Cell 3 analogue: load raw per-cell features from the cached parquet.

    The parquet concatenates all 185 per-subject SQLite joins. We filter to
    Erin's ``population`` (aggregated-CSV subjects minus the collapsed '370')
    and attach labels via a left join on ``subject``.
    """
    population = [s for s in ag_df["subject"].tolist() if s != "370"]
    cols = [
        "subject",
        "Nuclei_Location_Center_X",
        "Nuclei_Location_Center_Y",
        "Cells_AreaShape_Area",
        "Cells_AreaShape_MajorAxisLength",
        "Cells_AreaShape_MinorAxisLength",
        "Cells_AreaShape_Perimeter",
        "Nuclei_AreaShape_MajorAxisLength",
        "Nuclei_AreaShape_Area",
        "Cells_Intensity_MeanIntensity_Actin",
        "Nuclei_Intensity_MeanIntensity_Actin",
    ]
    df_new = pd.read_parquet(FIBROBLAST_SINGLECELL_PARQUET, columns=cols)
    df_new["subject"] = df_new["subject"].astype(str)
    df_new = df_new[df_new["subject"].isin(population)].copy()
    df_new = df_new.merge(ag_df[["subject", "label"]], on="subject", how="left")
    return df_new


def apply_erin_qc(df_new: pd.DataFrame) -> pd.DataFrame:
    """Cell 5 analogue: Erin's four QC gates in the original order."""
    print(f"Before border cell removal:        {df_new.shape}")
    mask_border = (
        (df_new["Nuclei_Location_Center_X"] > IM_WIDTH - BORDER_PX)
        | (df_new["Nuclei_Location_Center_X"] < BORDER_PX)
        | (df_new["Nuclei_Location_Center_Y"] > IM_HEIGHT - BORDER_PX)
        | (df_new["Nuclei_Location_Center_Y"] < BORDER_PX)
    )
    df = df_new.loc[~mask_border].reset_index(drop=True)
    print(f"After border cell removal:         {df.shape}")

    df["Cells2Nuclei_MajorAxisLengthRatio"] = df["Cells_AreaShape_MajorAxisLength"] / df["Nuclei_AreaShape_MajorAxisLength"]
    df["Cells2Nuclei_AreaShapeRatio"] = df["Cells_AreaShape_Area"] / df["Nuclei_AreaShape_Area"]
    df = df[df["Cells2Nuclei_MajorAxisLengthRatio"] > NUCLEI_MAJOR_RATIO_MIN].reset_index(drop=True)
    df = df[df["Cells2Nuclei_AreaShapeRatio"] > NUCLEI_AREA_RATIO_MIN].reset_index(drop=True)
    print(f"After cell seg prone cells removal:{df.shape}")

    df = df[
        (df["Cells_Intensity_MeanIntensity_Actin"] < CELLS_ACTIN_MAX)
        & (df["Nuclei_Intensity_MeanIntensity_Actin"] < NUCLEI_ACTIN_MAX)
    ].reset_index(drop=True)
    print(f"After intensity artifact removal:  {df.shape}")
    return df


def build_cellsize(df_1: pd.DataFrame) -> pd.DataFrame:
    """Cells 16+17 analogue: per-subject aggregation and micron conversion."""
    feat_cols_px = [
        "Cells_AreaShape_Area",
        "Cells_AreaShape_MajorAxisLength",
        "Cells_AreaShape_MinorAxisLength",
        "Cells_AreaShape_Perimeter",
    ]
    per_sub = df_1.groupby(["subject", "label"])[feat_cols_px].mean().reset_index()
    cellsize = per_sub.copy()
    cellsize["Area_Microns"] = cellsize["Cells_AreaShape_Area"] * (MICRONS_PER_PIXEL**2)
    cellsize["MajorAxisLength_Microns"] = cellsize["Cells_AreaShape_MajorAxisLength"] * MICRONS_PER_PIXEL
    cellsize["MinorAxisLength_Microns"] = cellsize["Cells_AreaShape_MinorAxisLength"] * MICRONS_PER_PIXEL
    cellsize["Perimeter_Microns"] = cellsize["Cells_AreaShape_Perimeter"] * MICRONS_PER_PIXEL

    psych = cellsize.loc[cellsize["label"].isin(["SZ", "SZA", "BP"])].copy()
    psych["label"] = "psychosis"
    return pd.concat([cellsize, psych], ignore_index=True)


def plot_boxplot(cellsize: pd.DataFrame, feat: str, out_path: Path) -> None:
    """Cell 18 analogue: one boxplot per feature, in Erin's exact style."""
    fig, axes = plt.subplots(1, 1, figsize=(6, 7))
    sns.boxplot(x="label", y=feat, data=cellsize, order=PATIENT_ORDER, ax=axes, palette=PALETTE, boxprops={"edgecolor": "none"})
    sns.swarmplot(x="label", y=feat, data=cellsize, order=PATIENT_ORDER, ax=axes, color=".2")

    mi, ma = cellsize[feat].describe()[["min", "max"]]
    for i, cat in enumerate(PATIENT_ORDER):
        if cat == "Control":
            continue
        pval = round_sig(
            ttest_ind(
                cellsize.loc[cellsize["label"] == "Control", feat].values,
                cellsize.loc[cellsize["label"] == cat, feat].values,
                equal_var=False,
            ).pvalue,
            2,
        )
        weight_f = "bold" if pval < 0.05 else None
        axes.annotate("p-value\n" + str(pval), xy=(i, mi - (0.07 * mi)), ha="center", size="small", weight=weight_f)

    axes.set_ylim([mi - (0.1 * mi), ma + (0.05 * ma)])
    axes.set(ylabel=feat + " (averaged per patient)")
    axes.set(xlabel="Patient category")
    plt.tight_layout()

    labels = [t.get_text() for t in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)

    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_pairplot(cellsize: pd.DataFrame, ag_df: pd.DataFrame, out_path: Path) -> None:
    """Cell 19 analogue: AreaShape-feature × slope pairplot.

    Erin computes slope from her recomputed df_mitoslope. Here we join on the
    ``slope`` column from the corrected aggregated CSV, which is the same
    slope the paper uses and is numerically aligned with Erin's recomputation
    after QC (the shipped CSV was produced by Marzieh's cell 4 from the same
    pipeline Erin reproduces).
    """
    slopes = ag_df[["subject", "label", "slope"]].copy()
    slopes["subject"] = slopes["subject"].astype(str)

    feat_cols_um = ["Area_Microns", "MajorAxisLength_Microns", "MinorAxisLength_Microns", "Perimeter_Microns"]
    corrplot = cellsize.merge(slopes, on=["subject", "label"], how="left")
    corrplot = corrplot.loc[corrplot["label"] != "psychosis"]
    pg = sns.pairplot(corrplot[feat_cols_um + ["slope"]])
    pg.savefig(out_path, bbox_inches="tight")
    plt.close(pg.fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    ag_df = load_corrected_ag_df()
    df_new = load_raw_single_cells(ag_df)
    df_1 = apply_erin_qc(df_new)
    cellsize = build_cellsize(df_1)

    feat_cols_um = [
        ("Area_Microns", "cell_area_microns.pdf"),
        ("MajorAxisLength_Microns", "major_axis_microns.pdf"),
        ("MinorAxisLength_Microns", "minor_axis_microns.pdf"),
        ("Perimeter_Microns", "perimeter_microns.pdf"),
    ]
    for feat, fname in feat_cols_um:
        plot_boxplot(cellsize, feat, OUT_DIR / fname)
        print(f"  wrote {OUT_DIR / fname}")

    plot_pairplot(cellsize, ag_df, OUT_DIR / "areashape_vs_slope_pairplot.pdf")
    print(f"  wrote {OUT_DIR / 'areashape_vs_slope_pairplot.pdf'}")

    print()
    print(f"Regenerated Erin's Supp Fig 3 (5 PDFs) from corrected labels → {OUT_DIR}")


if __name__ == "__main__":
    main()
