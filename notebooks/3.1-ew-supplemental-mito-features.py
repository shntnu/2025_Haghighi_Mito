# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: run_stuff
#     language: python
#     name: python3
# ---

# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

from haghighi_mito.config import AGGREGATED_PROFILES_PATH, PATIENT_ORDER, PATIENT_PALETTE, SUPPLEMENTAL_FIGURES_DIR


# %%
def round_sig(x, sig=2):
#     To round to n significant digits (non-zero digits)
    if x == 0:
        return 0
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


# %% [markdown]
# ## Read Data

# %%
df_1_avg_persub = pd.read_csv(AGGREGATED_PROFILES_PATH)

data_phs0 = df_1_avg_persub.groupby(['label','subject']).mean().reset_index()
data_phs_psych=data_phs0[data_phs0['label'].isin(["BP","SZ","SZA"])].copy()
data_phs_psych['label']='psychosis'
data_phs=pd.concat([data_phs0,data_phs_psych],ignore_index=True)

# %%
for targetFeature, nice_name, filename in [('Nuclei_ObjectSkeleton_NumberBranchEnds_mito_skel','Number of Branch Ends', 'branch_ends'),
                                 ('Nuclei_ObjectSkeleton_NumberBranchEnds_mito_skel', 'Number of Non Trunk Branches', 'non_trunk_branches'),
                                 ('Nuclei_ObjectSkeleton_NumberTrunks_mito_skel','Number of Trunks', 'num_trunks'),
                                 ('Nuclei_ObjectSkeleton_TotalObjectSkeletonLength_mito_skel','Total Mitochondria Skeleton Length', 'skeleton_length')]:
    fig, axes = plt.subplots(1, 1, figsize=(6,7))
    orderC = PATIENT_ORDER
    palette = PATIENT_PALETTE

    sns.boxplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes,palette=palette,boxprops={'edgecolor':'none'})
    sns.swarmplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes, color=".2")
    ############################
    mi,ma=data_phs[targetFeature].describe()[['min','max']]

    for i in range(len(orderC)):
        if orderC[i]!='Control':
            pval=round_sig(ttest_ind(data_phs.loc[data_phs['label']=='Control',targetFeature].values,\
                    data_phs.loc[data_phs['label']==orderC[i],targetFeature].values,equal_var=False).pvalue,2)
            print(pval)
            if pval<0.05:
                weight_f='bold'
            else:
                weight_f=None
            #axes.annotate('p-value\n'+str(pval),xy=(i,mi-0.12),horizontalalignment='center',size='small',weight=weight_f)
            axes.annotate('p-value\n'+str(pval),xy=(i,mi),horizontalalignment='center',size='small',weight=weight_f)

    axes.set(ylabel=f'{nice_name} (averaged per patient)')

    axes.set(xlabel='Patient category')
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)
    
    skel_dir = SUPPLEMENTAL_FIGURES_DIR / "skel_features"
    skel_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(skel_dir / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(skel_dir / f"{filename}.png", bbox_inches="tight")

# %%
for targetFeature, nice_name, filename in [('Cells_Intensity_IntegratedIntensity_Mito','Integrated Intensity of Mitotracker Fluorescence', 'integrated_intensity'),
                                 ('Cells_Intensity_MedianIntensity_Mito', 'Median Intensity of Mitotracker Fluorescence', 'median_intensity'),
                                 ('Cells_Intensity_MeanIntensity_Mito','Mean Intensity of Mitotracker Fluorescence', 'mean_intensity'),
                                 ('Cells_Intensity_StdIntensity_Mito','Standard Deviation of Intensity of Mitotracker Fluorescence', 'std_intensity')]:
    fig, axes = plt.subplots(1, 1, figsize=(6,7))
    orderC = PATIENT_ORDER
    palette = PATIENT_PALETTE

    sns.boxplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes,palette=palette,boxprops={'edgecolor':'none'})
    sns.swarmplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes, color=".2")
    ############################
    mi,ma=data_phs[targetFeature].describe()[['min','max']]

    for i in range(len(orderC)):
        if orderC[i]!='Control':
            pval=round_sig(ttest_ind(data_phs.loc[data_phs['label']=='Control',targetFeature].values,\
                    data_phs.loc[data_phs['label']==orderC[i],targetFeature].values,equal_var=False).pvalue,2)
            print(pval)
            if pval<0.05:
                weight_f='bold'
            else:
                weight_f=None
            axes.annotate('p-value\n'+str(pval),xy=(i,mi-0.25),horizontalalignment='center',size='small',weight=weight_f)
            #axes.annotate('p-value\n'+str(pval),xy=(i,mi),horizontalalignment='center',size='small',weight=weight_f)

    axes.set(ylabel=f'{nice_name} (averaged per patient)')
    axes.set_ylim([mi-0.5,ma+0.05])

    axes.set(xlabel='Patient category')
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)
    
    intensity_dir = SUPPLEMENTAL_FIGURES_DIR / "intensity_features"
    intensity_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(intensity_dir / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(intensity_dir / f"{filename}.png", bbox_inches="tight")

# %%
for targetFeature, nice_name, filename in [('Cells_Intensity_IntegratedIntensity_Mito','Integrated Intensity of Mitotracker Fluorescence in the Cell', 'integrated_intensity_cell'),
                                 ('Cells_Intensity_MedianIntensity_Mito', 'Median Intensity of Mitotracker Fluorescence in the Cell', 'median_intensity_cell'),
                                 ('Cells_Intensity_MeanIntensity_Mito','Mean Intensity of Mitotracker Fluorescence in the Cell', 'mean_intensity_cell'),
                                 ('Cells_Intensity_StdIntensity_Mito','Standard Deviation of Intensity of Mitotracker Fluorescence in the Cell', 'std_intensity_cell')]:
    fig, axes = plt.subplots(1, 1, figsize=(6,7))
    orderC = PATIENT_ORDER
    palette = PATIENT_PALETTE

    sns.boxplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes,palette=palette,boxprops={'edgecolor':'none'})
    sns.swarmplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes, color=".2")
    ############################
    mi,ma=data_phs[targetFeature].describe()[['min','max']]

    for i in range(len(orderC)):
        if orderC[i]!='Control':
            pval=round_sig(ttest_ind(data_phs.loc[data_phs['label']=='Control',targetFeature].values,\
                    data_phs.loc[data_phs['label']==orderC[i],targetFeature].values,equal_var=False).pvalue,2)
            print(pval)
            if pval<0.05:
                weight_f='bold'
            else:
                weight_f=None
            axes.annotate('p-value\n'+str(pval),xy=(i,mi-0.25),horizontalalignment='center',size='small',weight=weight_f)
            #axes.annotate('p-value\n'+str(pval),xy=(i,mi),horizontalalignment='center',size='small',weight=weight_f)

    axes.set(ylabel=f'{nice_name} (averaged per patient)')
    axes.set_ylim([mi-0.5,ma+0.05])

    axes.set(xlabel='Patient category')
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)
    
    intensity_dir = SUPPLEMENTAL_FIGURES_DIR / "intensity_features"
    intensity_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(intensity_dir / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(intensity_dir / f"{filename}.png", bbox_inches="tight")

# %%
for targetFeature, nice_name, filename in [('Cytoplasm_Intensity_IntegratedIntensity_Mito','Integrated Intensity of Mitotracker Fluorescence in the Cytoplasm', 'integrated_intensity_cytoplasm'),
                                 ('Cytoplasm_Intensity_MedianIntensity_Mito', 'Median Intensity of Mitotracker Fluorescence in the Cytoplasm', 'median_intensity_cytoplasm'),
                                 ('Cytoplasm_Intensity_MeanIntensity_Mito','Mean Intensity of Mitotracker Fluorescence in the Cytoplasm', 'mean_intensity_cytoplasm'),
                                 ('Cytoplasm_Intensity_StdIntensity_Mito','Standard Deviation of Intensity of Mitotracker Fluorescence in the Cytoplasm', 'std_intensity_cytoplasm')]:
    fig, axes = plt.subplots(1, 1, figsize=(6,7))
    orderC = PATIENT_ORDER
    palette = PATIENT_PALETTE

    sns.boxplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes,palette=palette,boxprops={'edgecolor':'none'})
    sns.swarmplot(x="label", y=targetFeature, data=data_phs,order=orderC,ax=axes, color=".2")
    ############################
    mi,ma=data_phs[targetFeature].describe()[['min','max']]

    for i in range(len(orderC)):
        if orderC[i]!='Control':
            pval=round_sig(ttest_ind(data_phs.loc[data_phs['label']=='Control',targetFeature].values,\
                    data_phs.loc[data_phs['label']==orderC[i],targetFeature].values,equal_var=False).pvalue,2)
            print(pval)
            if pval<0.05:
                weight_f='bold'
            else:
                weight_f=None
            axes.annotate('p-value\n'+str(pval),xy=(i,mi-0.25),horizontalalignment='center',size='small',weight=weight_f)
            #axes.annotate('p-value\n'+str(pval),xy=(i,mi),horizontalalignment='center',size='small',weight=weight_f)

    axes.set(ylabel=f'{nice_name} (averaged per patient)')
    axes.set_ylim([mi-0.5,ma+0.05])

    axes.set(xlabel='Patient category')
    plt.tight_layout()

    labels = [tick.get_text() for tick in axes.get_xticklabels()]
    labels[1] = "psychosis\n(BP + SZ + SZA)"
    axes.set_xticklabels(labels)
    
    intensity_dir = SUPPLEMENTAL_FIGURES_DIR / "intensity_features"
    intensity_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(intensity_dir / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(intensity_dir / f"{filename}.png", bbox_inches="tight")

# %%
