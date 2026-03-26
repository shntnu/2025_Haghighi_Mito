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
    PATIENT_FIGURES_DIR,
    PATIENT_LABELS_PATH,
    PATIENT_ORDER,
    PATIENT_PALETTE,
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
# ## Additional phenotype boxplots (Area, Perimeter, etc.)

# %%
plt.style.use(["science", "no-latex", "nature"])
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params, font_scale=1)

for tar in [
    "Cells_AreaShape_Area",
    "Cells_AreaShape_MajorAxisLength",
    "Cells_AreaShape_MinorAxisLength",
    "Cells_AreaShape_Perimeter",
    "Cells_MitoTubeness_MeanRadialPosition",
]:
    fig, axes = plt.subplots(1, 1, figsize=(6, 7))
    orderC = PATIENT_ORDER
    palette = PATIENT_PALETTE

    sns.boxplot(x="label", y=tar, data=data_phs, order=orderC, ax=axes, palette=palette, boxprops={"edgecolor": "none"})
    sns.swarmplot(x="label", y=tar, data=data_phs, order=orderC, ax=axes, color=".2")

    mi, ma = data_phs[tar].describe()[["min", "max"]]

    for i in range(len(orderC)):
        if orderC[i] != "Control":
            pval = round_sig(ttest_ind(data_phs.loc[data_phs["label"] == "Control", tar].values, data_phs.loc[data_phs["label"] == orderC[i], tar].values, equal_var=False).pvalue, 2)
            weight_f = "bold" if pval < 0.05 else None
            axes.annotate("p-value\n" + str(pval), xy=(i, mi - 0.12), horizontalalignment="center", size="small", weight=weight_f)

    axes.set_ylim([mi - 0.14, ma + 0.05])
    axes.set(ylabel=tar + " (averaged per patient)")
    axes.set(xlabel="Patient category")
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)

# %% [markdown]
# ## Choosing representative single cells for Figure 2

# %%
med = df_1_avg_persub.groupby("label")["slope"].transform("median")
df2 = df_1_avg_persub.assign(_dist=(df_1_avg_persub["slope"] - med).abs())

# One representative row per label: closest to the median slope
closest_to_median = df2.loc[df2.groupby("label")["_dist"].idxmin()].drop(columns=["_dist"]).sort_values("label")

closest_to_median[["subject", "label", "slope"]]
