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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from haghighi_mito.config import AGGREGATED_PROFILES_PATH, SUPPLEMENTAL_FIGURES_DIR

# %%
# total number of color bins to use for visualization
vis_bins = 20

nice_text = {'C':'Cell Measured with Bins Starting at the Center of the Cell',
             'CCN':'Cell Measured with Bins Starting at the Center of the Nucleus',
             'CEN':'Cell Measured with Bins Starting at the Edge of the Nucleus'}

# ordered to match other figures
custom_order = ['Control','BP','SZ','SZA','MDD']

# %%
df = pd.read_csv(AGGREGATED_PROFILES_PATH)

# filter to only mito_tubeness features
df = df[[x for x in df.columns if 'Cells_RadialDistribution' in x and 'mito_tubeness' in x] + ['label']]
df = df.groupby('label').median().reset_index()
df = df.set_index('label').reindex(custom_order)
df


# %% [markdown]
# # Functions

# %%
def rescale_df(df, vis_bins=False, castint=True):
    min_val = df.min().min()
    max_val = df.max().max()
    print(min_val, max_val)
    df = (df - min_val) / (max_val - min_val)
    min_val = df.min().min()
    max_val = df.max().max()
    print(min_val, max_val)
    if vis_bins:
        # scale from zero to vis_bins-1 for total bin number = vis_bins
        df = df*(vis_bins-1)
    if castint:
        df = df.astype(int)
    return df


# %% [markdown]
# # CEN, FracAtD, 12 bins - showing work

# %%
group = 'CEN'
nbins = 12

groupdf = df[[x for x in df.columns if f'of{nbins}' in x and f'_{group}_' in x]]
groupdf.columns = [x.replace(f'Cells_RadialDistribution_FracAtD_mito_tubeness_{group}_','').replace(f'of{nbins}','') for x in groupdf.columns]

# min-max rescale data, scale to number of visualization bins to work with visualization strategy
groupdf = rescale_df(groupdf, vis_bins=vis_bins, castint=True)
groupdf.reset_index(inplace=True)
groupdf

# %%
data_dict = {}
for label in groupdf['label']:
    data_dict[label] = [groupdf.loc[groupdf['label']==label,str(x)].item() for x in range(1,nbins+1)]
data_dict

# %%
# SET UP DATA FOR POLAR PLOTTING
# Create angles from 0 to 2*pi
theta = np.linspace(0, 2 * np.pi, 200).tolist()
ringintlist = []
ringlist = []
for bin in range(0, nbins):
    ringintlist += [bin]*len(theta)
    ringlist += [f'ring {bin}']*len(theta)
colorlist = sns.color_palette("YlOrBr", vis_bins)

palette_dict = {}
for label in data_dict.keys():
    palette ={}
    for bin in range(0, nbins):
        palette[f'ring {bin}'] = colorlist[data_dict[label][bin]]
    palette_dict[label] = palette

# %%
# make custom figure legend
min = df.min().min()
max = df.max().max()
mid = (max+min)/2
legend_elements = [Line2D([0], [0], color=colorlist[0], lw=4, label=f'{min:.2f}'),
                  Line2D([0], [0], color=colorlist[int(vis_bins/4)-1], lw=4, label=f'{(mid+min)/2:.2f}'),
                  Line2D([0], [0], color=colorlist[int(vis_bins/2)-1], lw=4, label=f'{mid:.2f}'),
                  Line2D([0], [0], color=colorlist[int(vis_bins*3/4)-1], lw=4, label=f'{(mid+max)/2:.2f}'),
                  Line2D([0], [0], color=colorlist[vis_bins-1], lw=4, label=f'{max:.2f}')]

# make subplots
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6*len(groupdf), 6),nrows=1, ncols=len(groupdf))

for i, label in enumerate(groupdf['label']):
    data = pd.DataFrame({'theta': theta*nbins, 'ringintlist': ringintlist,'ringlist': ringlist})
    g = sns.lineplot(data=data, x='theta', y='ringintlist', hue='ringlist', palette=palette_dict[label], ax=ax[i], linewidth=12)
    g.set(xlabel="")
    g.set(ylabel="")
    g.set(xticklabels=[])
    g.set(yticklabels=[])
    ax[i].set_title(label)
    ax[i].get_legend().remove()
fig.suptitle(f'{nbins} Bins, {nice_text[group]}', fontsize=16)
plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, .6))
# plt.savefig(f'/Users/eweisbar/Documents/projects/MitoDonna/paper/polar_plot_{nbins}_{group}.png', bbox_inches='tight')


# %% [markdown]
# # All conditions - mito tubeness features

# %%
for metric in ['FracAtD','MeanFrac','RadialCV']:
    metricdf = df[[x for x in df.columns if metric in x]]
    for group in ['CEN','CCN','C']:
        for nbins in [12,16]:
            groupdf = metricdf[[x for x in metricdf.columns if f'of{nbins}' in x and f'_{group}_' in x]]
            groupdf.columns = [x.replace(f'Cells_RadialDistribution_{metric}_mito_tubeness_{group}_','').replace(f'of{nbins}','') for x in groupdf.columns]

            # min-max rescale data, scale to number of visualization bins to work with visualization strategy
            groupdf = rescale_df(groupdf, vis_bins=vis_bins, castint=True)
            groupdf.reset_index(inplace=True)

            # Bin 1 is missing from CCN. Replace with zeros.
            if group =='CCN':
                groupdf['1'] = 0
            if group == 'CEN':
                # add 3 central rings to illustrate nucleus
                nbins += 3
                groupdf = groupdf.set_index('label')
                groupdf.columns = [str(int(x)+3) for x in groupdf.columns]
                groupdf = groupdf.reset_index()
                groupdf['1']  = groupdf['2'] = groupdf['3'] = 0

            data_dict = {}
            for label in groupdf['label']:
                data_dict[label] = [groupdf.loc[groupdf['label']==label,str(x)].item() for x in range(1,nbins+1)]
            
            # SET UP DATA FOR POLAR PLOTTING
            # Create angles from 0 to 2*pi
            theta = np.linspace(0, 2 * np.pi, 200).tolist()
            ringintlist = []
            ringlist = []
            for bin in range(0, nbins):
                ringintlist += [bin]*len(theta)
                ringlist += [f'ring {bin}']*len(theta)
            colorlist = sns.color_palette("YlOrBr", vis_bins)

            palette_dict = {}
            for label in data_dict.keys():
                palette ={}
                for bin in range(0, nbins):
                    palette[f'ring {bin}'] = colorlist[data_dict[label][bin]]
                if group == 'CEN':
                    # overwrite for 3 central rings added to signify nucleus
                    palette['ring 0'] = palette['ring 1'] = palette['ring 2'] = (0,0,1)    
                palette_dict[label] = palette
            

            # make custom figure legend
            min = df.min().min()
            max = df.max().max()
            mid = (max+min)/2
            legend_elements = [Line2D([0], [0], color=colorlist[0], lw=4, label=f'{min:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins/4)-1], lw=4, label=f'{(mid+min)/2:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins/2)-1], lw=4, label=f'{mid:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins*3/4)-1], lw=4, label=f'{(mid+max)/2:.2f}'),
                            Line2D([0], [0], color=colorlist[vis_bins-1], lw=4, label=f'{max:.2f}')]

            # make subplots
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6*len(groupdf), 6),nrows=1, ncols=len(groupdf))

            nrealbins = nbins
            if group == 'CEN':        
                #return nbins to number counted for labeling purposes
                nbins -= 3
            if nbins == 12:
                linewidth = 12.1
            elif nbins == 16:
                linewidth = 9
            for i, label in enumerate(groupdf['label']):
                data = pd.DataFrame({'theta': theta*nrealbins, 'ringintlist': ringintlist,'ringlist': ringlist})
                g = sns.lineplot(data=data, x='theta', y='ringintlist', hue='ringlist', palette=palette_dict[label], ax=ax[i], linewidth=linewidth)
                if group == 'CEN':
                    # fill small empty spot in the middle
                    sns.scatterplot(data=[0,0], color='blue', s=100, ax=ax[i])
                else:
                    sns.scatterplot(data=[0,0], color=palette_dict[label]['ring 0'], s=150, ax=ax[i])
                g.set(xlabel="",ylabel="",xticklabels=[],yticklabels=[])
                ax[i].set_title(label)
                ax[i].get_legend().remove()
                ax[i].grid(False)
                fig.suptitle(f'{nbins} Bins, {nice_text[group]}', fontsize=16)
            plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, .6),title=metric)
            skel_metric_dir = SUPPLEMENTAL_FIGURES_DIR / "skeleton" / metric
            skel_metric_dir.mkdir(parents=True, exist_ok=True)
            plt.savefig(skel_metric_dir / f"polar_plot_{nbins}_{group}.png", bbox_inches='tight')

# %% [markdown]
# # Intensity features in Mito channel

# %%
df = pd.read_csv(AGGREGATED_PROFILES_PATH)
# filter to only mito_tubeness features
df = df[[x for x in df.columns if 'Mito' in x and 'Radial' in x and 'Nuclei' not in x] + ['label']]
df = df.groupby('label').median().reset_index()
df = df.set_index('label').reindex(custom_order)
df


# %%
nrealbins = 4
for metric in ['FracAtD','MeanFrac','RadialCV']:
    metricdf = df[[x for x in df.columns if metric in x]]
    for group, ndrawbins in [('Cytoplasm',7), ('Cells',4)]:
            groupdf = metricdf[[x for x in metricdf.columns if f'of{nrealbins}' in x and f'{group}' in x]]
            groupdf.columns = [x.replace(f'{group}_RadialDistribution_{metric}_Mito_','').replace(f'of{nrealbins}','') for x in groupdf.columns]

            # min-max rescale data, scale to number of visualization bins to work with visualization strategy
            groupdf = rescale_df(groupdf, vis_bins=vis_bins, castint=True)
            groupdf.reset_index(inplace=True)

            if group == 'Cytoplasm':
                groupdf = groupdf.set_index('label')
                groupdf.columns = [str(int(x)+3) for x in groupdf.columns]
                groupdf = groupdf.reset_index()
                groupdf['1']  = groupdf['2'] = groupdf['3'] = 0

            data_dict = {}
            for label in groupdf['label']:
                data_dict[label] = [groupdf.loc[groupdf['label']==label,str(x)].item() for x in range(1,ndrawbins+1)]
            
            # SET UP DATA FOR POLAR PLOTTING
            # Create angles from 0 to 2*pi
            theta = np.linspace(0, 2 * np.pi, 200).tolist()
            ringintlist = []
            ringlist = []
            for bin in range(0, ndrawbins):
                ringintlist += [bin]*len(theta)
                ringlist += [f'ring {bin}']*len(theta)
            colorlist = sns.color_palette("YlOrBr", vis_bins)

            palette_dict = {}
            for label in data_dict.keys():
                palette ={}
                for bin in range(0, ndrawbins):
                    palette[f'ring {bin}'] = colorlist[data_dict[label][bin]]
                if group == 'Cytoplasm':
                    # overwrite for 3 central rings added to signify nucleus
                    palette['ring 0'] = palette['ring 1'] = palette['ring 2'] = (0,0,1)    
                palette_dict[label] = palette
            

            # make custom figure legend
            min = df.min().min()
            max = df.max().max()
            mid = (max+min)/2
            legend_elements = [Line2D([0], [0], color=colorlist[0], lw=4, label=f'{min:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins/4)-1], lw=4, label=f'{(mid+min)/2:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins/2)-1], lw=4, label=f'{mid:.2f}'),
                            Line2D([0], [0], color=colorlist[int(vis_bins*3/4)-1], lw=4, label=f'{(mid+max)/2:.2f}'),
                            Line2D([0], [0], color=colorlist[vis_bins-1], lw=4, label=f'{max:.2f}')]

            # make subplots
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6*len(groupdf), 6),nrows=1, ncols=len(groupdf))

            if group == 'Cells':
                linewidth = 50
            elif group == 'Cytoplasm':
                linewidth = 25
            for i, label in enumerate(groupdf['label']):
                data = pd.DataFrame({'theta': theta*ndrawbins, 'ringintlist': ringintlist,'ringlist': ringlist})
                g = sns.lineplot(data=data, x='theta', y='ringintlist', hue='ringlist', palette=palette_dict[label], ax=ax[i], linewidth=linewidth)
                # fill empty spot in the middle
                if group == 'Cytoplasm':
                    sns.scatterplot(data=[0,0], color='blue', s=1000, ax=ax[i])
                elif group == 'Cells':
                    sns.scatterplot(data=[0,0], color=palette_dict[label]['ring 0'], s=2200, ax=ax[i])
                g.set(xlabel="",ylabel="",xticklabels=[],yticklabels=[])
                ax[i].set_title(label)
                ax[i].get_legend().remove()
                ax[i].grid(False)
                fig.suptitle(f'{nrealbins} Bins, {group}', fontsize=16)
            plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, .6),title=metric)
            fluor_metric_dir = SUPPLEMENTAL_FIGURES_DIR / "fluorescence" / metric
            fluor_metric_dir.mkdir(parents=True, exist_ok=True)
            plt.savefig(fluor_metric_dir / f"polar_plot_{nrealbins}_{group}.png", bbox_inches='tight')

# %%
