# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Patient MITO-SLOPE Analysis
#
# Ported from upstream `phenotype_discovery/2.slope_analysis.ipynb`.
#
# This notebook generates manuscript figures for patient fibroblast mitochondrial
# dispersion analysis (Figure 4b, Figure 4c) and performs feature correlation analysis.
#
# **Data source:** Uses pre-aggregated patient profiles (`aggregated_profiles_fibroblast.csv`)
# downloaded from S3. The heavy single-cell processing pipeline that produces this CSV
# (SQLite → QC → standardize → aggregate) is documented but not run here.

# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
import scipy.stats as stats
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import ttest_ind

from haghighi_mito.config import (
    AGGREGATED_PROFILES_PATH,
    FIBROBLAST_DATA_DIR,
    PATIENT_FIGURES_DIR,
    PATIENT_LABELS_PATH,
    PATIENT_ORDER,
    PATIENT_PALETTE,
    PIXEL_SIZE_UM,
    SUPPLEMENTAL_FIGURES_DIR,
)

# Ensure output directories exist
PATIENT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# %%
def round_sig(x, sig=2):
    """Round to n significant digits."""
    if x == 0:
        return 0
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


# %% [markdown]
# ## Load Pre-Aggregated Patient Data

# %%
df_1_avg_persub = pd.read_csv(AGGREGATED_PROFILES_PATH)

# Apply subject label corrections (from upstream notebook)
df_1_avg_persub.loc[df_1_avg_persub["subject"].astype(str) == "272", "label"] = "Control"
df_1_avg_persub.loc[df_1_avg_persub["subject"].astype(str) == "MCL004", "label"] = "SZA"
df_1_avg_persub["label"] = df_1_avg_persub["label"].replace("MDD or Dep", "MDD")

print(f"Loaded {len(df_1_avg_persub)} patient profiles")
print(f"Labels: {df_1_avg_persub['label'].value_counts().to_dict()}")

# %%
# Prepare psychosis composite group (BP + SZ + SZA)
data_phs0 = df_1_avg_persub.groupby(["label", "subject"]).mean(numeric_only=True).reset_index()
data_phs_psych = data_phs0[data_phs0["label"].isin(["BP", "SZ", "SZA"])].copy()
data_phs_psych["label"] = "psychosis"
data_phs = pd.concat([data_phs0, data_phs_psych], ignore_index=True)

# %% [markdown]
# ## Figure 4c — MITO-SLOPE per patient category

# %%
plt.style.use(["science", "no-latex", "nature"])
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params, font_scale=1)

targetFeature = "slope"

fig, axes = plt.subplots(1, 1, figsize=(6, 7))
orderC = PATIENT_ORDER
palette = PATIENT_PALETTE

sns.boxplot(x="label", y=targetFeature, data=data_phs, order=orderC, ax=axes, palette=palette, boxprops={"edgecolor": "none"})
sns.swarmplot(x="label", y=targetFeature, data=data_phs, order=orderC, ax=axes, color=".2")

mi, ma = data_phs[targetFeature].describe()[["min", "max"]]

for i in range(len(orderC)):
    if orderC[i] != "Control":
        pval = round_sig(ttest_ind(data_phs.loc[data_phs["label"] == "Control", targetFeature].values, data_phs.loc[data_phs["label"] == orderC[i], targetFeature].values, equal_var=False).pvalue, 2)
        print(pval)
        weight_f = "bold" if pval < 0.05 else None
        axes.annotate("p-value\n" + str(pval), xy=(i, mi - 0.12), horizontalalignment="center", size="small", weight=weight_f)

axes.set_ylim([mi - 0.14, ma + 0.05])

# Bracket annotations for psychosis group
ofset1 = 1
ofset2 = 1.04
axes.annotate("", xy=(3, 1.199 - ofset1), xytext=(3, 1.20 - ofset1), fontsize=4, arrowprops=dict(arrowstyle="-[, widthB=20.0, lengthB=2", color="#916e99", lw=2.0))
axes.annotate("", xy=(2, 1.259 - ofset2), xytext=(2, 1.26 - ofset2), fontsize=4, arrowprops=dict(arrowstyle="-[, widthB=18, lengthB=3", color="#916e99", lw=2.0))
axes.text(x=0.5, y=0.18, s="psychosis", color="#916e99", fontsize=11)
axes.text(x=1.6, y=0.16, s="BP        +        SZ       +      SZA", color="#916e99", fontsize=11)

axes.set(ylabel="Mitochondrial dispersion (averaged per patient)")
axes.set(xlabel="Patient category")
plt.tight_layout()

labels = [tick.get_text() for tick in axes.get_xticklabels()]
labels[1] = "psychosis\n(BP + SZ + SZA)"
axes.set_xticklabels(labels)

fig.savefig(PATIENT_FIGURES_DIR / "Figure4c.pdf", bbox_inches="tight")
print(f"Saved {PATIENT_FIGURES_DIR / 'Figure4c.pdf'}")

# %% [markdown]
# ## Figure 4b — Radial distribution profiles by patient category
#
# This figure requires the per-cell standardized data (not just aggregated profiles).
# The pre-aggregated CSV contains per-subject means, but Figure 4b needs per-cell
# distributions for error bars. Skipping for now — requires single-cell data.

# %% [markdown]
# ## Weighted Mitochondrial Radial Position

# %%
data_phs["Cells_MitoTubeness_MeanRadialPosition"] = 0.0

for x in range(1, 13):
    data_phs["Cells_MitoTubeness_MeanRadialPosition"] += (x / 12) * data_phs[f"Cells_RadialDistribution_MeanFrac_mito_tubeness_CCN_{x}of12"]

# %% [markdown]
# ## MeanRadialPosition boxplot (z-scored, unit-less diagnostic)
#
# MeanRadialPosition is a fraction (0–1) weighted-sum of radial bins, so no
# micron conversion applies. Uses the z-scored aggregated CSV.

# %%
plt.style.use(["science", "no-latex", "nature"])
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params, font_scale=1)

SUPPLEMENTAL_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

tar = "Cells_MitoTubeness_MeanRadialPosition"
fig, axes = plt.subplots(1, 1, figsize=(6, 7))
sns.boxplot(x="label", y=tar, data=data_phs, order=PATIENT_ORDER, ax=axes, palette=PATIENT_PALETTE, boxprops={"edgecolor": "none"})
sns.swarmplot(x="label", y=tar, data=data_phs, order=PATIENT_ORDER, ax=axes, color=".2")

mi, ma = data_phs[tar].describe()[["min", "max"]]

for i in range(len(PATIENT_ORDER)):
    if PATIENT_ORDER[i] != "Control":
        pval = round_sig(ttest_ind(data_phs.loc[data_phs["label"] == "Control", tar].values, data_phs.loc[data_phs["label"] == PATIENT_ORDER[i], tar].values, equal_var=False).pvalue, 2)
        weight_f = "bold" if pval < 0.05 else None
        axes.annotate("p-value\n" + str(pval), xy=(i, mi - 0.12), horizontalalignment="center", size="small", weight=weight_f)

axes.set_ylim([mi - 0.14, ma + 0.05])
axes.set(ylabel=tar + " (averaged per patient)")
axes.set(xlabel="Patient category")
plt.tight_layout()

labels = [tick.get_text() for tick in axes.get_xticklabels()]
labels[1] = "psychosis\n(BP + SZ + SZA)"
axes.set_xticklabels(labels)

fig.savefig(SUPPLEMENTAL_FIGURES_DIR / f"patient_boxplot_{tar}_zscored.pdf", bbox_inches="tight")
plt.close(fig)

# %% [markdown]
# ## Supp Fig 3E — Feature list and correlation values
# Correlations are computed on the raw-micron pair plot below (panel E).

# %%
panel_e_features = [
    "Cells_AreaShape_Area",
    "Cells_AreaShape_MajorAxisLength",
    "Cells_AreaShape_MinorAxisLength",
    "Cells_AreaShape_Perimeter",
]

# %% [markdown]
# ## Supp Fig 3 A–D — Raw-micron AreaShape boxplots
#
# Loads per-cell raw pixel values from the allFeatures pickle, applies subject
# label corrections, aggregates per subject, converts pixels → microns using
# ``PIXEL_SIZE_UM`` from config (currently provisional — see note there).

# %%
pkl_path = FIBROBLAST_DATA_DIR / "single_cell_with_annot_allFeatures.pkl"
df_cells = pd.read_pickle(pkl_path, compression="infer")
df_cells["subject"] = df_cells["subject"].replace(["370E", "370F", "370H"], "370")
df_cells.loc[df_cells["subject"].astype(str) == "272", "label"] = "Control"
df_cells.loc[df_cells["subject"].astype(str) == "MCL004", "label"] = "SZA"
df_cells["label"] = df_cells["label"].replace("MDD or Dep", "MDD")

raw_area_features = [
    "Cells_AreaShape_Area",
    "Cells_AreaShape_MajorAxisLength",
    "Cells_AreaShape_MinorAxisLength",
    "Cells_AreaShape_Perimeter",
]
data_raw = df_cells.groupby(["label", "subject"])[raw_area_features].mean().reset_index()
data_raw_psych = data_raw[data_raw["label"].isin(["BP", "SZ", "SZA"])].copy()
data_raw_psych["label"] = "psychosis"
data_raw = pd.concat([data_raw, data_raw_psych], ignore_index=True)

# Pixel → micron conversion (area uses squared factor)
for feat in ["Cells_AreaShape_MajorAxisLength", "Cells_AreaShape_MinorAxisLength", "Cells_AreaShape_Perimeter"]:
    data_raw[feat] = data_raw[feat] * PIXEL_SIZE_UM
data_raw["Cells_AreaShape_Area"] = data_raw["Cells_AreaShape_Area"] * (PIXEL_SIZE_UM**2)

# (feature_key, nice_name, unit, short_slug) — matches ErinWeisbart's
# plot_mito_features.ipynb convention (nice names + short filenames).
panel_ad_features = [
    ("Cells_AreaShape_Area", "Cell area", "µm²", "cell_area"),
    ("Cells_AreaShape_Perimeter", "Cell perimeter", "µm", "cell_perimeter"),
    ("Cells_AreaShape_MajorAxisLength", "Major axis length", "µm", "major_axis"),
    ("Cells_AreaShape_MinorAxisLength", "Minor axis length", "µm", "minor_axis"),
]

plt.style.use(["science", "no-latex", "nature"])
sns.set_theme(style="ticks", rc=custom_params, font_scale=1)


for feat, nice_name, unit, slug in panel_ad_features:
    fig, axes = plt.subplots(1, 1, figsize=(6, 7))
    sns.boxplot(x="label", y=feat, data=data_raw, order=PATIENT_ORDER, ax=axes, palette=PATIENT_PALETTE, boxprops={"edgecolor": "none"})
    sns.swarmplot(x="label", y=feat, data=data_raw, order=PATIENT_ORDER, ax=axes, color=".2")

    mi, ma = data_raw[feat].describe()[["min", "max"]]
    for i in range(len(PATIENT_ORDER)):
        if PATIENT_ORDER[i] != "Control":
            pval = round_sig(ttest_ind(data_raw.loc[data_raw["label"] == "Control", feat].values, data_raw.loc[data_raw["label"] == PATIENT_ORDER[i], feat].values, equal_var=False).pvalue, 2)
            weight_f = "bold" if pval < 0.05 else None
            axes.annotate("p-value\n" + str(pval), xy=(i, mi), ha="center", size="small", weight=weight_f)

    axes.set(ylabel=f"{nice_name} ({unit}, averaged per patient)")
    axes.set(xlabel="Patient category")
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)

    fig.savefig(SUPPLEMENTAL_FIGURES_DIR / f"SuppFig3_{slug}_microns.pdf", bbox_inches="tight")
    plt.close(fig)
print(f"Saved raw-micron SuppFig3 A–D to {SUPPLEMENTAL_FIGURES_DIR} (PIXEL_SIZE_UM={PIXEL_SIZE_UM})")

# %% [markdown]
# ## Supp Fig 3E — Raw-micron pair plot
#
# Uses raw-micron AreaShape values aggregated from the pickle, joined with the
# z-scored ``slope`` from the aggregated CSV (slope is derived from radial
# features, not AreaShape, so leaving it z-scored is fine).

# %%
raw_persub = df_cells.groupby("subject")[raw_area_features].mean()
for feat in ["Cells_AreaShape_MajorAxisLength", "Cells_AreaShape_MinorAxisLength", "Cells_AreaShape_Perimeter"]:
    raw_persub[feat] = raw_persub[feat] * PIXEL_SIZE_UM
raw_persub["Cells_AreaShape_Area"] = raw_persub["Cells_AreaShape_Area"] * (PIXEL_SIZE_UM**2)

slope_persub = df_1_avg_persub.set_index("subject")["slope"]
panel_e_raw = raw_persub.join(slope_persub, how="inner").reset_index()

panel_e_raw_corr = panel_e_raw[panel_e_features].corrwith(panel_e_raw["slope"])
print("Raw-micron AreaShape correlations with MITO-SLOPE:")
print(panel_e_raw_corr)

pg = sns.pairplot(panel_e_raw[panel_e_features + ["slope"]])
for ax in pg.axes.flat:
    if ax is not None:
        ax.set_xlabel(ax.get_xlabel(), fontsize=8)
        ax.set_ylabel(ax.get_ylabel(), fontsize=8)
        ax.tick_params(axis="both", labelsize=6)
corr_text = "\n".join(f"{feat.split('_')[-1]}: corr with MITO-SLOPE = {corr_val:.3f}" for feat, corr_val in panel_e_raw_corr.items())
pg.figure.suptitle(corr_text, x=0.5, y=1.02, ha="center", va="bottom", fontsize=10)
pg.figure.savefig(SUPPLEMENTAL_FIGURES_DIR / "SuppFig3E_AreaShape_vs_slope_pairplot_microns.pdf", bbox_inches="tight")
plt.close(pg.figure)
print(f"Saved raw-micron SuppFig3E pairplot (n={len(panel_e_raw)})")

# %% [markdown]
# ## Supp Fig 3 A–E — Combined figure (gdraw layout)
#
# Single-PDF assembly matching the current Google Drawing:
#   - Top: 2×2 grid of boxplots (A=Area, B=Perimeter, C=MajorAxis, D=MinorAxis)
#   - Bottom: 5×5 pairplot (E) with correlation text column on the right
# Panel labels A–E are drawn in bold. Category axis shows n (patients / group).

# %%
panel_boxplots = [
    ("A", "Cells_AreaShape_Area", "Cell area", "µm²"),
    ("B", "Cells_AreaShape_Perimeter", "Cell perimeter", "µm"),
    ("C", "Cells_AreaShape_MajorAxisLength", "Major axis length", "µm"),
    ("D", "Cells_AreaShape_MinorAxisLength", "Minor axis length", "µm"),
]

# n = unique patients per group (for x-tick annotation)
n_per_group = data_raw.groupby("label")["subject"].nunique().to_dict()
xlabels_with_n = [f"{g}\nn = {n_per_group.get(g, 0)}" if g != "psychosis" else f"psychosis\n(BP+SZ+SZA)\nn = {n_per_group.get(g, 0)}" for g in PATIENT_ORDER]

plt.style.use(["science", "no-latex", "nature"])
sns.set_theme(style="ticks", rc=custom_params, font_scale=1)

fig = plt.figure(figsize=(14, 16))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.3], hspace=0.28)

# Top: 2×2 boxplot grid
top = outer[0].subgridspec(2, 2, hspace=0.45, wspace=0.25)
for idx, (letter, feat, short_title, unit) in enumerate(panel_boxplots):
    ax = fig.add_subplot(top[idx // 2, idx % 2])
    sns.boxplot(x="label", y=feat, data=data_raw, order=PATIENT_ORDER, ax=ax, palette=PATIENT_PALETTE, boxprops={"edgecolor": "none"})
    sns.swarmplot(x="label", y=feat, data=data_raw, order=PATIENT_ORDER, ax=ax, color=".2", size=3)

    mi, ma = data_raw[feat].describe()[["min", "max"]]
    for i, grp in enumerate(PATIENT_ORDER):
        if grp != "Control":
            pval = round_sig(ttest_ind(data_raw.loc[data_raw["label"] == "Control", feat].values, data_raw.loc[data_raw["label"] == grp, feat].values, equal_var=False).pvalue, 2)
            weight_f = "bold" if pval < 0.05 else None
            ax.annotate(f"p={pval}", xy=(i, mi), ha="center", size=8, weight=weight_f)

    ax.set_ylabel(f"{short_title} ({unit})", fontsize=10)
    ax.set_xlabel("")
    ax.set_xticklabels(xlabels_with_n, fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    # Panel letter
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")

# Bottom: pairplot (5x5) + correlation text column
bottom = outer[1].subgridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
pair_grid = bottom[0].subgridspec(len(panel_e_features) + 1, len(panel_e_features) + 1, hspace=0.1, wspace=0.1)
cols_e = panel_e_features + ["slope"]
short_names = {
    "Cells_AreaShape_Area": "Area",
    "Cells_AreaShape_MajorAxisLength": "Major axis",
    "Cells_AreaShape_MinorAxisLength": "Minor axis",
    "Cells_AreaShape_Perimeter": "Perimeter",
    "slope": "MITO-SLOPE",
}
n_e = len(cols_e)
for i, yf in enumerate(cols_e):
    for j, xf in enumerate(cols_e):
        ax = fig.add_subplot(pair_grid[i, j])
        if i == j:
            ax.hist(panel_e_raw[xf].dropna(), bins=20, color="#4c72b0", edgecolor="white", linewidth=0.3)
        else:
            ax.scatter(panel_e_raw[xf], panel_e_raw[yf], s=8, color="#4c72b0", alpha=0.7, edgecolor="none")
        ax.tick_params(axis="both", labelsize=6)
        if i < n_e - 1:
            ax.set_xticklabels([])
        if j > 0:
            ax.set_yticklabels([])
        if i == n_e - 1:
            ax.set_xlabel(short_names[xf], fontsize=8)
        if j == 0:
            ax.set_ylabel(short_names[yf], fontsize=8)

# Panel E letter (outside top-left of pairplot)
ax_e_letter = fig.add_subplot(bottom[0])
ax_e_letter.set_axis_off()
ax_e_letter.text(-0.08, 1.02, "E", transform=ax_e_letter.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")

# Correlation text column on right
ax_corr = fig.add_subplot(bottom[1])
ax_corr.set_axis_off()
lines = [f"{short_names[f]}:\ncorrelation with\nMITO-SLOPE\n= {panel_e_raw_corr[f]:.3f}\n" for f in panel_e_features]
ax_corr.text(0.02, 0.95, "\n".join(lines), transform=ax_corr.transAxes, fontsize=10, va="top", ha="left")

combined_path = SUPPLEMENTAL_FIGURES_DIR / "SuppFig3_combined_microns.pdf"
fig.savefig(combined_path, bbox_inches="tight")
plt.close(fig)
print(f"Saved combined figure: {combined_path}")
