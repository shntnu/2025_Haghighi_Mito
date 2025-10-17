# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: gan
#     language: python
#     name: gan
# ---

# %% [markdown]
# ## Pipeline to screen phenotype strength of a target feature in various datasets
# * this pipeline is for screening per site values of the target feature versus control wells for a given dataset

# %%
# %load_ext autoreload
# %autoreload 2
# %matplotlib notebook
import os
import pickle
import time
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import sklearn.preprocessing as sp

today = date.today()

import sys

from singlecell.preprocess import (
    extract_cpfeature_names,
    find_highly_correlated_features,
    handle_nans,
)
from singlecell.preprocess.control_for_cellcount import control_feature_y_for_variable_x
from singlecell.preprocess.filter_out_edge_single_cells import edgeCellFilter
from singlecell.process import bbf_test, normalize_funcs, precision_recall, statistical_tests

# singlecell-morph is now installed via uv, no need for sys.path.insert
from singlecell.read import read_single_cell_sql
from singlecell.save.save_pandas_dfs import saveDF_to_CSV_GZ_no_timestamp
from singlecell.visualize import visualize_n_SingleCell

# %% [markdown]
#
#
#
# # Steps:
#
# ### 1. Form a list of orthogonal features to target feature based on per well aggregated profiles
#    - Read per well csv profiles and form a df for all plates
#    - The selected set of features have <0.1 correlation with the target feature according to the target dataset.
#    - Save the list as a pickle file
#      - Output Folder on s3: <font color='blue'>projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/results/target_orth_features_lists</font>
#
# ### 2. Create per site aggregate level of data for a few features
#    - Read saved orthogonal to target feature list which is derived from step 1
#    - Read the set of orth+target features from each plate and form per site measures
#    - Concat and save it in output folder
#      - Output Folder on s3: <font color='blue'>projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/per_site_aggregated_profiles</font>
#
# ### 3. Load per_site aggregated data and control target feature for cell counts
#   - Read the saved aggregated per site level data
#   - PER PLATE control of target feature for cell count
#   - PER PLATE low variance feature removal
#
#
# ### 4. Apply PER PLATE statistical test between the target feature values for each pert versus controls
# - And the same for orth features
# - We calculate test stats per plate and take the p and t values of the plate with median t value
#     - Output Folder on s3: <font color='blue'>projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/results/reverse_phenotype_strength</font>
#

# %%
############################################################################################
# This list contains features which their phenotype strength were significantly different across
# various categories in the patient fibroblast data and were available in at least one of the
# three datasets examined here

# target_features_list_orf_cdrp=['Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16',\
#              'Nuclei_Texture_DifferenceVariance_Mito_10_00_256',\
#              'Nuclei_Texture_Contrast_Mito_10_00_256']

target_features_list_orf_cdrp = ["slope"]
target_features_list_lincs = ["slope"]
# target_features_list_lincs=['Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16']

########################## Project root directory and path to results ########################
# home_path="/home/ubuntu/" # ec2
home_path = "/home/jupyter-mhaghigh@broadinst-ee45a/"  # dgx
mito_project_root_dir = (
    home_path + "bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/"
)
save_results_dir = mito_project_root_dir + "/workspace/results/"

# %%
# python3 ~/imaging-backup-scripts/restore_intelligent.py

# %%
# https://imaging-platform.s3.amazonaws.com/projects/2018_04_20_Rosetta/workspace/preprocessed_data

# %% [markdown]
# #### Preprocess datasets metadata into a unified format to read

# %%
# ########## jump_orf
dataset = "jump_orf"
jump_orf_meta_tsv = mito_project_root_dir + "workspace/metadata/JUMP-ORF/ORF_list.tsv"
annot = pd.read_csv(jump_orf_meta_tsv, sep="\t")
annot["batch_plate"] = annot["Batch"] + "-" + annot["Assay_Plate_Barcode"]
annot["ctrl_well"] = annot["Symbol"].isin(["LacZ", "BFP", "HcRed", "LUCIFERASE"])
# annot['ctrl_well']=annot['Symbol'].isin(['LacZ'])
# annot['pert_id']=annot['broad_sample']
annot["Metadata_pert_type"] = annot["pert_type"]
annot.to_csv(mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv")

#
# # ########## lincs
dataset = "lincs"
# annot=pd.read_csv("/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/raw-profiles/CP_LINCS/metadata/matadata_lincs.csv")
annot0 = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/lincs/DrugRepurposing_Metadata.csv"
)
annot = pd.read_csv(mito_project_root_dir + "/workspace/metadata/LINCS_meta.csv")
annot = annot.merge(
    annot0[["Metadata_Plate", "Metadata_Well", "Metadata_pert_name"]],
    how="left",
    on=["Metadata_Plate", "Metadata_Well"],
)
annot["Batch"] = "2016_04_01_a549_48hr_batch1_Mito_Project"
annot["batch_plate"] = annot["Batch"] + "-" + annot["Metadata_Plate"]
annot["ctrl_well"] = annot["Metadata_pert_type"].isin(["control"])
# annot['pert_id']=annot['Metadata_pert_id_dose']
annot.to_csv(mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv")

# # ########### CDRP
# dataset="CDRP"
# # annot=pd.read_csv("/home/ubuntu/gallery/cpg0012-wawer-bioactivecompoundprofiling/broad/workspace/metadata/platemaps/CDRP/barcode_platemap.csv")
# annot=pd.read_csv(mito_project_root_dir+"/workspace/metadata/CDRP_meta.csv")
# annot['Batch']='CDRP'
# annot['batch_plate']=annot['Batch']+'-'+annot['Metadata_Plate'].astype(str)
# annot['ctrl_well']=annot['Metadata_pert_type'].isin(['control'])
# # annot['pert_id']=annot['broad_sample']
# annot.to_csv(mito_project_root_dir+"/workspace/metadata/preprocessed/annot_"+dataset+'.csv',index=False)

# ########## jump_orf/jump_crispr/jump_compound

plates = pd.read_csv(mito_project_root_dir + "/workspace/metadata/JUMP/plate.csv.gz")
wells = pd.read_csv(mito_project_root_dir + "/workspace/metadata/JUMP/well.csv.gz")
compound = pd.read_csv(mito_project_root_dir + "/workspace/metadata/JUMP/compound.csv.gz")
orf = pd.read_csv(mito_project_root_dir + "/workspace/metadata/JUMP/orf.csv.gz")
crispr = pd.read_csv(mito_project_root_dir + "/workspace/metadata/JUMP/crispr.csv.gz")

compound_plates = plates[plates["Metadata_PlateType"] == "COMPOUND"].reset_index(drop=True)

dataset = "jump_orf"
annot_orf = wells.merge(orf, on=["Metadata_JCP2022"]).merge(
    plates, on=["Metadata_Plate", "Metadata_Source"]
)
annot_orf["Batch"] = annot_orf["Metadata_Batch"]
annot_orf["batch_plate"] = annot_orf["Metadata_Batch"] + "-" + annot_orf["Metadata_Plate"]
annot_orf["ctrl_well"] = annot_orf["Metadata_Symbol"].isin(["LacZ", "BFP", "HcRed", "LUCIFERASE"])
annot_orf.to_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv"
)

dataset = "jump_crispr"
annot_crispr = wells.merge(crispr, on=["Metadata_JCP2022"]).merge(
    plates, on=["Metadata_Plate", "Metadata_Source"]
)
annot_compound = wells.merge(compound, on=["Metadata_JCP2022"]).merge(
    compound_plates, on=["Metadata_Plate", "Metadata_Source"]
)

annot_crispr["Batch"] = annot_crispr["Metadata_Batch"]
annot_crispr["batch_plate"] = annot_crispr["Metadata_Batch"] + "-" + annot_crispr["Metadata_Plate"]
annot_compound["batch_plate"] = (
    annot_compound["Metadata_Batch"] + "-" + annot_compound["Metadata_Plate"]
)

annot_crispr["ctrl_well"] = annot_crispr["Metadata_Symbol"].isin(["non-targeting"])
annot_crispr.to_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv"
)

dataset = "jump_compound"
# annot['ctrl_well']=annot['Symbol'].isin(['LacZ'])
annot_compound.to_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv"
)


# ########### TA-ORF
dataset = "taorf"
# annot=pd.read_csv("/home/ubuntu/gallery/cpg0012-wawer-bioactivecompoundprofiling/broad/workspace/metadata/platemaps/CDRP/barcode_platemap.csv")
annot_taorf = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz"
)

annot_taorf["Batch"] = "2013_10_11_SIGMA2_Pilot"
annot_taorf["batch_plate"] = annot_taorf["Batch"] + "-" + annot_taorf["Metadata_Plate"].astype(str)
annot_taorf["ctrl_well"] = annot_taorf["Metadata_gene_name"].isin(["LacZ", "Luciferase"])
annot_taorf["Metadata_pert_type"] = annot_taorf["Metadata_gene_name"].isin(["LacZ", "Luciferase"])

## annot_taorf.Metadata_ASSAY_WELL_ROLE.unique()  array(['Untreated', 'Treated', 'CTRL'], dtype=object)
annot_taorf["Metadata_pert_type"] = annot_taorf["Metadata_ASSAY_WELL_ROLE"]

annot_taorf[
    [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_gene_name",
        "Metadata_pert_name",
        "Metadata_pert_type",
        "Metadata_broad_sample",
        "Metadata_moa",
        "batch_plate",
        "Batch",
        "ctrl_well",
    ]
].to_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv",
    index=False,
)

# %%
# https://cellpainting-gallery.s3.amazonaws.com/cpg0003-rosetta/broad/workspace/preprocessed_data

# %%
# annot_taorf=pd.read_csv(mito_project_root_dir+"/workspace/metadata/TA-ORF/replicate_level_cp_normalized.csv.gz")
annot_taorf = pd.read_csv(
    "/home/ubuntu/bucket/projects/2018_04_20_Rosetta/workspace/preprocessed_data/TA-ORF-BBBC037-Rohban/CellPainting/replicate_level_cp_augmented.csv.gz"
)

# %%
# aws s3 sync s3://imaging-platform/projects/2018_04_20_Rosetta/workspace/preprocessed_data/ s3://imaging-platform/projects/2018_04_20_Rosetta/workspace/preprocessed_data_gallery/

# %%
# # ls -R /home/ubuntu/gallery/cpg0003-rosetta/broad/workspace/preprocessed_data

# %%
# https://imaging-platform.s3.amazonaws.com/projects/2018_04_20_Rosetta/workspace/preprocessed_data/TA-ORF-BBBC037-Rohban/CellPainting/replicate_level_cp_normalized.csv

# %%
# annot_taorf['Metadata_broad_sample'].unique()

# %%
# Metadata_broad_sample

# %%
# annot_crispr[annot_crispr["Metadata_Symbol"].str.contains('GSK3')]

# %%
annot_taorf[annot_taorf["Metadata_broad_sample"] == "DMSO"].Metadata_ASSAY_WELL_ROLE.unique()

# %%
annot_taorf.Metadata_ASSAY_WELL_ROLE.unique()

# %%
# annot_taorf[annot_taorf["Metadata_gene_name"].isin(["LacZ","Luciferase"])]
# annot_taorf[annot_taorf.columns[annot_taorf.columns.str.contains("Metadata")]]
# [['Metadata_Plate','Metadata_Well','Metadata_gene_name','Metadata_pert_name','Metadata_broad_sample','Metadata_moa']]

# %%
# ls /home/ubuntu/gallery/cpg0017-rohban-pathways/broad/workspace/backend/2013_10_11_SIGMA2_Pilot/

# %%

# %% [markdown]
# ### Set dataset specific parameters

# %%
lincs_meta_cols = [
    "Metadata_broad_sample",
    "Metadata_dose_recode",
    "Metadata_pert_id",
    "Metadata_pert_mfc_id",
    "Metadata_InChIKey14",
    "Metadata_pert_type",
    "Metadata_moa",
    "Metadata_target",
    "Metadata_pert_id_dose",
    "Metadata_pert_name",
]

# lincs_meta_cols=['Metadata_broad_sample','Metadata_dose_recode','Metadata_pert_id','Metadata_pert_mfc_id',\
# 'Metadata_InChIKey14','Metadata_pert_type','Metadata_pert_id_dose']

cdrp_meta_cols = [
    "Metadata_broad_sample",
    "Metadata_mmoles_per_liter2",
    "Metadata_pert_id",
    "Metadata_Sample_Dose",
    "Metadata_moa",
]
jumporf_meta_cols = ["Metadata_Symbol", "Metadata_broad_sample", "Metadata_JCP2022"]
jumpcrispr_meta_cols = ["Metadata_NCBI_Gene_ID", "Metadata_Symbol", "Metadata_JCP2022"]
jumpcompound_meta_cols = ["Metadata_InChIKey", "Metadata_InChI", "Metadata_JCP2022"]
taorf_meta_cols = [
    "Metadata_gene_name",
    "Metadata_pert_name",
    "Metadata_broad_sample",
    "Metadata_moa",
]


# jump_orf_params={'profiles_path':"/home/ubuntu/jumpbucket/projects/2021_04_26_Production/workspace/backend/",\
#                  'meta_cols':jumporf_meta_cols,\
#                  'pert_col':'broad_sample',\
#                  'target_features_list':target_features_list_orf_cdrp
#                 }

jump_orf_params = {
    "profiles_path": home_path + "/gallery/cpg0016-jump/source_4/workspace/backend/",
    "meta_cols": jumporf_meta_cols,
    "pert_col": "Metadata_JCP2022",
    "target_features_list": target_features_list_orf_cdrp,
}

cdrp_params = {
    "profiles_path": home_path
    + "/gallery/cpg0012-wawer-bioactivecompoundprofiling/broad/workspace/backend/",
    "meta_cols": cdrp_meta_cols,
    "pert_col": "Metadata_Sample_Dose",
    "target_features_list": target_features_list_orf_cdrp,
}

lincs_params = {
    "profiles_path": home_path
    + "/bucket/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace/backend/",
    "meta_cols": lincs_meta_cols,
    "pert_col": "Metadata_pert_id_dose",
    "target_features_list": target_features_list_lincs,
}

jump_crispr_params = {
    "profiles_path": home_path + "/gallery/cpg0016-jump/source_13/workspace/backend/",
    "meta_cols": jumpcrispr_meta_cols,
    "pert_col": "Metadata_JCP2022",
    "target_features_list": target_features_list_orf_cdrp,
}

jump_compound_params = {
    "profiles_path": home_path + "/gallery/cpg0016-jump/source_13/workspace/backend/",
    "meta_cols": jumpcompound_meta_cols,
    "pert_col": "Metadata_JCP2022",
    "target_features_list": target_features_list_orf_cdrp,
}


ta_orf_params = {
    "profiles_path": "~/gallery/cpg0017-rohban-pathways/broad/workspace/backend/",
    "meta_cols": taorf_meta_cols,
    "pert_col": "Metadata_broad_sample",
    "target_features_list": target_features_list_orf_cdrp,
}

ds_info_dict = {
    "jump_orf": jump_orf_params,
    "CDRP": cdrp_params,
    "lincs": lincs_params,
    "jump_crispr": jump_crispr_params,
    "jump_compound": jump_compound_params,
    "taorf": ta_orf_params,
}


# 'broad_sample', 'pert_type', 'control_type'

# results=annot[['Symbol','broad_sample', 'pert_type', 'control_type']].drop_duplicates().reset_index(drop=True)

# dataset='CDRP';dataset_meta_hue='Metadata_moa'
# dataset='lincs';dataset_meta_hue='Metadata_moa'
# dataset='jump_orf';dataset_meta_hue='Symbol'

# %%

# %%
import gc
from functools import reduce

import pandas as pd
from sqlalchemy import create_engine


def read_per_well_data(
    input_data_dir,
    annot,
    prof_workspace_folder_name="profiles",
    fformat=".parquet",
):
    batches = annot["Batch"].unique()

    df_agg_all_batches_ls = []
    for b in batches:
        print(b)
        #         if "Metadata_Source" in annot.columns:
        source_str = annot.loc[annot["Batch"] == b, "Metadata_Source"].unique()[0]
        #             print(source_str)
        profile_path = (
            input_data_dir + source_str + "/workspace/" + prof_workspace_folder_name + "/"
        )
        #         else:
        #             profile_path = input_data_dir + "/workspace/profiles/"

        df_sag_ls = []
        plates_exist = os.listdir(profile_path + b)
        plates_meta = annot.loc[annot["Batch"] == b, "Metadata_Plate"].unique()
        plates = set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)

            fileName = profile_path + b + "/" + p + "/" + p + fformat
            #             print(fileName)
            if os.path.exists(fileName):
                if fformat == ".parquet":
                    sc_df = pd.read_parquet(fileName)
                elif fformat in [".csv", ".csv.gz"]:
                    sc_df = pd.read_csv(fileName)

                #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
                sc_df["Metadata_Batch"] = b
                sc_df["Metadata_Plate"] = p
                df_sag_ls.append(sc_df)
                del sc_df
                gc.collect()
            else:
                print(fileName, " not exists")

        if df_sag_ls:
            df_sag = pd.concat(df_sag_ls, axis=0)
            df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches = pd.concat(df_agg_all_batches_ls, axis=0, ignore_index=True)
    return df_agg_all_batches


def read_per_well_data_csvs(input_data_dir, annot):
    batches = annot["Batch"].unique()

    df_agg_all_batches_ls = []
    for b in batches:
        print(b)
        df_sag_ls = []
        plates_exist = os.listdir(input_data_dir + b)
        plates_meta = annot.loc[annot["Batch"] == b, "Metadata_Plate"].unique()
        plates = set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)

            fileName = input_data_dir + b + "/" + p + "/" + p + ".csv"
            #             print(fileName)
            if os.path.exists(fileName):
                sc_df = pd.read_csv(fileName)

                #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
                sc_df["Metadata_Batch"] = b
                sc_df["Metadata_Plate"] = p
                df_sag_ls.append(sc_df)
                del sc_df
                gc.collect()
            else:
                print(fileName, " not exists")

        if df_sag_ls:
            df_sag = pd.concat(df_sag_ls, axis=0)
            df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches = pd.concat(df_agg_all_batches_ls, axis=0)
    return df_agg_all_batches


def sample_single_cells_from_sql(input_data_dir, annot):
    batches = annot["Batch"].unique()

    df_agg_all_batches_ls = []
    for b in batches:
        print(b)
        df_sag_ls = []
        plates_exist = os.listdir(input_data_dir + b)
        plates_meta = annot.loc[annot["Batch"] == b, "Metadata_Plate"].unique()
        plates = list(set(plates_meta) & set(plates_exist))
        for p in plates[:20]:
            fileName = input_data_dir + b + "/" + p + "/" + p + ".sqlite"
            print(p, fileName)
            n_rand_ims = 100
            sc_df = read_single_cell_sql.readSingleCellData_sqlalch_random_image_subset(
                fileName, n_rand_ims
            )
            #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
            sc_df["Metadata_Batch"] = b
            sc_df["Metadata_Plate"] = p
            df_sag_ls.append(sc_df)
            del sc_df
            gc.collect()

        df_sag = pd.concat(df_sag_ls, axis=0)
        df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches = pd.concat(df_agg_all_batches_ls, axis=0)
    return df_agg_all_batches


def form_per_site_aggregated_profiles(annot, input_data_dir, output_dir, feature_list2):
    import gc

    import pandas as pd

    batches = annot["Batch"].unique().tolist()[113:]
    for b in batches:
        df_sag_ls = []
        src = annot.loc[annot["Batch"] == b, "Metadata_Source"].unique()[0]
        print(b, src)
        input_data_dir = "/".join([src if "source_" in i else i for i in input_data_dir.split("/")])

        plates_exist = os.listdir(input_data_dir + b)
        plates_meta = annot.loc[annot["Batch"] == b, "Metadata_Plate"].unique()
        plates = set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)
            fileName = input_data_dir + b + "/" + p + "/" + p + ".sqlite"
            print(fileName)
            sc_df = read_single_cell_sql.readSingleCellData_sqlalch_features_subset(
                fileName, feature_list2
            )
            cell_count_col_name = sc_df.columns[sc_df.columns.str.contains("Count_Cell")].values[0]
            per_site_aggregate = (
                sc_df.groupby(["Metadata_Well", "Metadata_Site"])
                .mean(numeric_only=True)[feature_list2 + [cell_count_col_name]]
                .reset_index()
            )
            #                                                                         ['Count_Cells']].reset_index()
            #                                                                   ['Count_CellsIncludingEdges']].reset_index() for crispr

            per_site_aggregate["Count_Cells"] = per_site_aggregate[cell_count_col_name]
            per_site_aggregate["Metadata_Batch"] = b
            per_site_aggregate["Metadata_Plate"] = p
            df_sag_ls.append(per_site_aggregate)
        #             del sc_df
        #             gc.collect()

        df_sag = pd.concat(df_sag_ls, axis=0)
        fileNameToSave = output_dir + "/" + b + "_site_agg_profiles"
        print(fileNameToSave)
        saveDF_to_CSV_GZ_no_timestamp(df_sag, fileNameToSave)

    return


from scipy.stats import f


def TwoSampleT2Test(X, Y):
    nx, p = X.shape
    ny, _ = Y.shape
    delta = np.mean(X, axis=0) - np.mean(Y, axis=0)
    Sx = np.cov(X, rowvar=False)
    Sy = np.cov(Y, rowvar=False)
    S_pooled = ((nx - 1) * Sx + (ny - 1) * Sy) / (nx + ny - 2)
    S_pooled = S_pooled + np.eye(S_pooled.shape[0]) * 1e-6
    t_squared = (
        (nx * ny)
        / (nx + ny)
        * np.matmul(np.matmul(delta.transpose(), np.linalg.inv(S_pooled)), delta)
    )
    statistic = t_squared * (nx + ny - p - 1) / (p * (nx + ny - 2))
    F = f(p, nx + ny - p - 1)
    p_value = 1 - F.cdf(statistic)
    #     print(f"Test statistic: {statistic}\nDegrees of freedom: {p} and {nx+ny-p-1}\np-value: {p_value}")

    # Convert F-statistic to z-score
    z_score = (statistic - (p / (nx + ny - p - 1))) / np.sqrt(
        (2 * p * (nx + ny - p - 1)) / ((nx + ny - 2) * (nx + ny - p - 1))
    )
    std_p_val = 2 * (1 - norm.cdf(abs(z_score)))

    return statistic, p_value, std_p_val


from scipy.stats import chi2


def HotellingsT_internal(X, Y, test="f"):
    n1, p = X.shape
    n2 = Y.shape[0]

    mu = np.zeros(p)

    # Calculate means and differences
    Xmeans = np.mean(X, axis=0)
    Ymeans = np.mean(Y, axis=0)
    X_diff = X - Xmeans
    Y_diff = Y - Ymeans

    # Calculate pooled covariance matrix
    S_pooled = 1 / (n1 + n2 - 2) * (X_diff.T @ X_diff + Y_diff.T @ Y_diff)

    # Calculate test statistic
    diff_means = Xmeans - Ymeans - mu
    if test == "f":
        test_statistic = (
            n1
            * n2
            / (n1 + n2)
            * diff_means
            @ np.linalg.inv(S_pooled)
            @ diff_means.T
            * (n1 + n2 - p - 1)
            / (p * (n1 + n2 - 2))
        )
        df1 = p
        df2 = n1 + n2 - p - 1
        p_value = 1 - f.cdf(test_statistic, df1, df2)
    elif test == "chi":
        test_statistic = n1 * n2 / (n1 + n2) * diff_means @ np.linalg.inv(S_pooled) @ diff_means.T
        df1 = p
        df2 = None
        p_value = 1 - chi2.cdf(test_statistic, df1)
    else:
        return "Invalid test type"

    return test_statistic, p_value


from scipy.signal import find_peaks


# NOTE: Renamed from find_end_slope to avoid duplicate definition at line 673.
# This earlier version uses width=2 for peak detection and returns (np.nan, np.nan) when no extrema found,
# while the active version (line 673) uses width=1 and returns (np.nan, 0).
def find_end_slope_width2(data, height=None):
    peaks, _ = find_peaks(data, height=height, width=2)
    valleys, _ = find_peaks(-data, height=height, width=2)
    extermas = np.concatenate((peaks, valleys))
    if extermas.size == 0:
        return np.nan, np.nan

    last_peak_ind = np.max(extermas)
    slope = data[-1] - data[last_peak_ind]
    return last_peak_ind, slope


# def find_end_slope2(data, height=None):
#     min_max_indc = [np.argmax(data), np.argmin(data)]
#     last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]
#     if last_peak_ind0 == []:
#         return 0, 0
#     last_peak_ind = np.max(last_peak_ind0)
#     slope = (data[-1] - data[last_peak_ind]) / (len(data) - last_peak_ind)

#     #     last_peak_ind = 0
#     #     slope = data[-1] - data[-4]
#     return last_peak_ind, slope
# # def find_peaks_valleys(data, height=None):
# #     peaks, _ = find_peaks(data, height=height,width=2)
# #     valleys, _ = find_peaks(-data, height=height,width=2)
# #     return peaks, valleys


from scipy.signal import savgol_filter


def smooth_data(data, window_length=5, polyorder=3):
    return savgol_filter(data, window_length, polyorder)


def find_end_slope(data, height=None):
    peaks, _ = find_peaks(data, height=height, width=1)
    valleys, _ = find_peaks(-data, height=height, width=1)
    extermas = np.concatenate((peaks, valleys))
    if extermas.size == 0:
        return np.nan, 0

    last_peak_ind = np.max(extermas)
    slope = data[-1] - data[last_peak_ind]
    return last_peak_ind, slope


def subtract_control(group):
    batch_plate = group.name
    control_values = control_df_perplate.loc[batch_plate]
    return group - control_values


def find_end_slope2(data, height=None):
    data = smooth_data(data)
    #     min_max_indc = [np.argmax(data[3:] + 3), np.argmin(data[3:] + 3)]
    min_max_indc = [np.argmax(data), np.argmin(data)]
    last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]
    if last_peak_ind0 == []:
        return 0, 0
    last_peak_ind = np.max(last_peak_ind0)
    last_two_points_amplitude = (data[-1] + data[-2]) / 2
    slope = (last_two_points_amplitude - data[last_peak_ind]) / (len(data) - last_peak_ind - 1)

    #     last_peak_ind = 0
    #     slope = data[-2]  # - data[-2]
    return last_peak_ind, slope


# def find_end_slope2(data, height=None, plot=False):
#     min_max_indc = [np.argmax(data), np.argmin(data)]
#     last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]
#     if last_peak_ind0 == []:
#         return 0, 0
#     last_peak_ind = np.max(last_peak_ind0)
#     slope = (data[-1] - data[last_peak_ind]) / (len(data) - last_peak_ind-1)

#     # If plot flag is enabled, create plot
#     if plot:
#         x_values = range(len(data))
#         plt.plot(x_values, data, label="Data", color="blue")
#         y_values_slope = [data[last_peak_ind] + slope * (x - last_peak_ind) for x in x_values]
#         plt.plot(x_values, y_values_slope, label="Slope", color="red")
#         plt.legend()
#         plt.show()

#     return last_peak_ind, slope


# %%
dataset = "lincs"
dataset = "jump_orf"

##################### Read preprocessed metadata
annot = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv",
    dtype={"Metadata_Plate": str},
)

annot

# %%
annot

# %%
annot[(annot["Metadata_Symbol"].isnull())].Metadata_broad_sample.unique()

# %%
# annot#["pert_type"].unique()

# %%
# df_sag.columns[df_sag.columns.str.contains('Cyto.*MeanFrac')]

# %% [markdown] heading_collapsed=true
# ## 1- form a list of orthogonal features based on per well aggregated profiles
# - jump-orf: 40  orth features saved!
# - CDRP:     23  orth features saved!
# - lincs:    161  orth features (threshold=0.9), 96  orth features (threshold=0.8) saved!
# - jump_crispr:    14  orth features saved!
# - taorf:    14  orth features saved!

# %% hidden=true
# %time

# dataset='jump_orf'
# dataset='CDRP'
dataset = "lincs"
# dataset="jump_crispr"
# dataset="taorf"

f_substr = "MeanFrac"
##################### Read preprocessed metadata
annot = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv",
    dtype={"Metadata_Plate": str},
)
# target_features_list=ds_info_dict[dataset]["target_features_list"]+['slope']
# target_features_list=

##################### Read per well data
if dataset == "lincs":  # for this batch of lincs data we dont have well level profiles
    df_agg_all_batches = sample_single_cells_from_sql(ds_info_dict[dataset]["profiles_path"], annot)
else:
    df_agg_all_batches = read_per_well_data_csvs(ds_info_dict[dataset]["profiles_path"], annot)
    df_agg_all_batches = read_per_well_data(
        ds_info_dict[dataset]["profiles_path"],
        # annot_source, # FIXME: annot_source is not defined
        annot,
        ds_info_dict[dataset]["prof_workspace_folder_name"],
        fformat=ds_info_dict[dataset]["pformat"],
    )


# ##################### Clean and shrink features
cp_features, cp_features_analysis_0 = extract_cpfeature_names.extract_cpfeature_names(
    df_agg_all_batches
)
df_sag, cp_features_analysis = handle_nans.handle_nans(
    df_agg_all_batches,
    cp_features_analysis_0,
    thrsh_null_ratio=0.05,
    thrsh_std=0.001,
    fill_na_method="drop-rows",
)
# df_sag['batch_plate']=df_sag['Metadata_Batch']+df_sag['Metadata_Plate']


##################################### merge with annot
common_cols_2merge = list(set(annot.columns) & set(df_sag.columns))
df_sag["Metadata_Plate"] = df_sag["Metadata_Plate"].astype(str)
df_sag = pd.merge(df_sag, annot, how="left", on=common_cols_2merge)
df_sag = df_sag[~df_sag["batch_plate"].isnull()].reset_index(drop=True)

df_sag = normalize_funcs.standardize_per_catX(df_sag, "batch_plate", cp_features_analysis).copy()

df_sag = df_sag[~df_sag["ctrl_well"].isnull()].reset_index(drop=True)


##########################
f_substr = "MeanFrac"
# f_substr='RadialCV'

target_columns = [
    "Cells_RadialDistribution_" + f_substr + "_mito_tubeness_" + str(i) + "of16"
    for i in range(5, 17)
]
df_ctrl_targetFs = (
    df_sag[df_sag["ctrl_well"]]
    .groupby("batch_plate")[target_columns]
    .mean()
    .reset_index(drop=False)
)
df_sag_2 = pd.merge(df_sag, df_ctrl_targetFs, how="left", on="batch_plate")

diff_pattern_arr = (
    df_sag_2[[t + "_x" for t in target_columns]].values
    - df_sag_2[[t + "_y" for t in target_columns]].values
)
x = np.apply_along_axis(find_end_slope2, 1, diff_pattern_arr)

df_sag_2["slope"] = x[:, 1]
##################### Calculate correlation of target feature with the rest of features

# target_features_list=[t+'_x' for t in target_columns]+
target_features_list = [
    "slope",
    "Cytoplasm_RadialDistribution_MeanFrac_mito_tubeness_16of16",
    "Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16",
]


uncorr_with_tfs = set()
for tfeat in target_features_list:
    if tfeat in df_sag_2.columns:
        corr_math = (
            df_sag_2[list(set(cp_features_analysis) & set(df_sag_2.columns))]
            .corrwith(df_sag_2[tfeat], axis=0)
            .abs()
            .reset_index()
        )
        if len(uncorr_with_tfs) != 0:
            uncorr_with_tfs = set(uncorr_with_tfs) & set(
                corr_math.loc[corr_math[0] < 0.1, "index"].tolist()
            )
        else:
            uncorr_with_tfs = set(corr_math.loc[corr_math[0] < 0.1, "index"].tolist())


uncorr_with_tfs = [
    t for t in uncorr_with_tfs if "RadialDistribution" not in t and "_Correlation" not in t
]
# uncorr_with_tfs=[t for t in uncorr_with_tfs if "RadialDistribution" not in t]
##################### remove correlated features in the orth feature list to remove redundancy
# and increase efficency in computations
similar_fs_2remove = find_highly_correlated_features.find_correlation(
    df_sag_2[uncorr_with_tfs], threshold=0.6, remove_negative=True
)
uncorr_feats_condese = list(set(uncorr_with_tfs) - set(similar_fs_2remove))


pd.DataFrame({"orth_fs": uncorr_feats_condese}).to_csv(
    save_results_dir + "target_pattern_orth_features_lists/" + dataset + "_2.csv", index=False
)
print(len(uncorr_feats_condese), " orth features saved!")

# %% hidden=true
# df_agg_all_batches=sample_single_cells_from_sql(ds_info_dict[dataset]["profiles_path"],annot);

# %% hidden=true
# df_agg_all_batches

# %% hidden=true
pd.DataFrame({"orth_fs": uncorr_feats_condese}).to_csv(
    save_results_dir + "target_pattern_orth_features_lists/" + dataset + "_2.csv", index=False
)
print(len(uncorr_feats_condese), " orth features saved!")


# %% hidden=true

# %% hidden=true
# len(uncorr_with_tfs)
# # corr_math
# uncorr_feats_condese
# df_sag[df_sag['batch_plate'].isnull()]

# %% hidden=true
similar_fs_2remove = find_highly_correlated_features.find_correlation(
    df_sag_2[uncorr_with_tfs], threshold=0.8, remove_negative=True
)
uncorr_feats_condese = list(set(uncorr_with_tfs) - set(similar_fs_2remove))
pd.DataFrame({"orth_fs": uncorr_feats_condese}).to_csv(
    save_results_dir + "target_pattern_orth_features_lists/" + dataset + ".csv", index=False
)
print(len(uncorr_feats_condese), " orth features saved!")

# %% hidden=true
# annot[]

# %% hidden=true
print(len(uncorr_feats_condese), " orth features saved!"), dataset

# %% hidden=true
# tfeat in df_sag_2.columns

# %% hidden=true
target_features_list = ds_info_dict[dataset]["target_features_list"] + ["slope"]

# %% hidden=true
# df_sag_2
target_features_list

# %% hidden=true
# b.shape,diff_pattern_arr.shape,x.shape

# %% hidden=true
# output_save_dir
save_results_dir

# %% [markdown] heading_collapsed=true
# ## 2. Create per site aggregate level of data for a few features
# - Fix a subset of features
#   - Target feature + its orthogonal features
#
# - Read the subset of features from each plate and form per site measures
# - Concat and save it in output folder

# %% hidden=true
dataset = "CDRP"
# dataset='lincs'
# dataset="taorf"
dataset = "jump_orf"
# dataset="jump_crispr"
dataset = "jump_compound"  # source_2 excluded
##################### Read preprocessed metadata
annot = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv",
    dtype={"Metadata_Plate": str},
)
# target_features_list=ds_info_dict[dataset]["target_features_list"]

output_save_dir = mito_project_root_dir + "/workspace/per_site_aggregated_profiles_newpattern_2/"

if 1:
    if 0:
        uncorr_feats_condese_1 = pd.read_csv(
            save_results_dir + "target_pattern_orth_features_lists/fibroblast_derived.csv"
        )["orth_fs"].tolist()
        uncorr_feats_condese_2 = pd.read_csv(
            save_results_dir + "target_pattern_orth_features_lists/" + dataset + ".csv"
        )["orth_fs"].tolist()
        uncorr_feats_condese = list(set(uncorr_feats_condese_1 + uncorr_feats_condese_2))
    else:
        if dataset == "lincs":
            uncorr_feats_condese = [
                "Nuclei_AreaShape_FormFactor",
                "Nuclei_AreaShape_Eccentricity",
                "Cells_AreaShape_Solidity",
                "Cells_Intensity_MaxIntensity_Mito",
                "Cells_AreaShape_Eccentricity",
                "Cytoplasm_AreaShape_MaxFeretDiameter",
                "Nuclei_Texture_AngularSecondMoment_DNA_8_45",
            ]
        else:
            uncorr_feats_condese = pd.read_csv(
                save_results_dir + "target_pattern_orth_features_lists/fibroblast_derived.csv"
            )["orth_fs"].tolist()

else:
    uncorr_feats_condese = pd.read_csv(
        save_results_dir + "target_pattern_orth_features_lists/" + dataset + ".csv"
    )["orth_fs"].tolist()

radial_meanFrac_features = [
    "Cells_RadialDistribution_MeanFrac_mito_tubeness_" + str(i) + "of16" for i in range(1, 17)
]

feature_list2 = uncorr_feats_condese + radial_meanFrac_features

output_dir = output_save_dir + dataset
form_per_site_aggregated_profiles(
    annot, ds_info_dict[dataset]["profiles_path"], output_dir, feature_list2
)
# %%
uncorr_feats_condese = pd.read_csv(
    save_results_dir + "target_pattern_orth_features_lists/fibroblast_derived.csv"
)["orth_fs"].tolist()

# %% hidden=true
len(uncorr_feats_condese)

# %% hidden=true
# ls /home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/per_site_aggregated_profiles/jump_crispr

# %% [markdown] heading_collapsed=true
# ## 3. Load per_site aggregated data and control target feature for cell counts
#   - Read the saved aggregated per site level data
#   - PER PLATE control of target feature for cell count
#   - PER PLATE low variance feature removal

# %% hidden=true
import pandas as pd
from sklearn import linear_model

f_substr = "MeanFrac"
target_columns = [
    "Cells_RadialDistribution_" + f_substr + "_mito_tubeness_" + str(i) + "of16"
    for i in range(5, 17)
]

dataset = "CDRP"
dataset = "jump_orf"
# dataset="lincs"
dataset = "jump_crispr"
dataset = "taorf"
dataset = "jump_compound"

per_site_profiles_path = (
    mito_project_root_dir + "/workspace/per_site_aggregated_profiles_newpattern_2/"
)
# feature_list=['Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16',\
#              'Nuclei_Texture_DifferenceVariance_Mito_10_00_256',\
#              'Nuclei_Texture_Contrast_Mito_10_00_256']
# feature_list=['Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16']

if 0:
    uncorr_feats_condese = pd.read_csv(
        save_results_dir + "target_pattern_orth_features_lists/" + dataset + ".csv"
    )["orth_fs"].tolist()

else:
    if dataset == "lincs":
        uncorr_feats_condese = [
            "Nuclei_AreaShape_FormFactor",
            "Nuclei_AreaShape_Eccentricity",
            "Cells_AreaShape_Solidity",
            "Cells_Intensity_MaxIntensity_Mito",
            "Cells_AreaShape_Eccentricity",
            "Cytoplasm_AreaShape_MaxFeretDiameter",
            "Nuclei_Texture_AngularSecondMoment_DNA_8_45",
        ]
    else:
        uncorr_feats_condese = pd.read_csv(
            save_results_dir + "target_pattern_orth_features_lists/fibroblast_derived.csv"
        )["orth_fs"].tolist()


annot = pd.read_csv(
    mito_project_root_dir + "/workspace/metadata/preprocessed/annot_" + dataset + ".csv",
    dtype={"Metadata_Plate": str},
)

target_features_list = ds_info_dict[dataset]["target_features_list"]

cols2remove_lowVars_eachPlate = []
per_site_df_ls = []

batches = annot["Batch"].unique()  # [:5]
# print(annot['Batch'].unique().shape)
# sdfdsfds
for b in batches:
    fileNameToSave = (
        per_site_profiles_path + "/" + dataset + "/" + b + "_site_agg_profiles" + ".csv.gz"
    )

    if os.path.exists(fileNameToSave):
        per_site_df_b = pd.read_csv(fileNameToSave)

        #     for tfeat in target_features_list:
        #         per_site_df_b=control_feature_y_for_variable_x(per_site_df_b,tfeat,'Count_Cells','_ccOut')

        thrsh_std = 0.001
        cols2remove_lowVars_eachPlate += (
            per_site_df_b[uncorr_feats_condese]
            .std()[per_site_df_b[uncorr_feats_condese].std() < thrsh_std]
            .index.tolist()
        )

        per_site_df_ls.append(per_site_df_b)

per_site_df = pd.concat(per_site_df_ls, axis=0, ignore_index=True)

per_site_df, cp_features_analysiss = handle_nans.handle_nans(
    per_site_df,
    target_columns + uncorr_feats_condese,
    thrsh_null_ratio=0.05,
    thrsh_std=0.001,
    fill_na_method="drop-rows",
)
# per_site_df['batch_plate']=per_site_df['Metadata_Batch']+'-'+per_site_df['Metadata_Plate']


# per_site_df=pd.merge(per_site_df, annot, how='left',left_on=['batch_plate','Metadata_Well'],\
#                      right_on=['batch_plate','Metadata_Well'])


common_cols_2merge = list(set(annot.columns) & set(per_site_df.columns))
# annot['Metadata_Plate'] = annot['Metadata_Plate'].astype(str)
per_site_df["Metadata_Plate"] = per_site_df["Metadata_Plate"].astype(str)
# per_site_df['Metadata_Batch'] = per_site_df['Metadata_Batch'].astype(str)
# per_site_df['Metadata_Well'] = per_site_df['Metadata_Well'].astype(str)

merge_how = "inner" if dataset == "jump_crispr" or dataset == "jump_compound" else "left"

per_site_df = pd.merge(per_site_df, annot, how=merge_how, on=common_cols_2merge)

if "Metadata_pert_type" in per_site_df.columns:
    per_site_df = per_site_df[~per_site_df["Metadata_pert_type"].isnull()].reset_index(drop=True)

# uncorr_feats_condese = per_site_df.columns[per_site_df.columns.str.contains("Cells_|Nuclei_|Cytoplasm_")].tolist()[1:]
uncorr_feats_cond = list(set(uncorr_feats_condese) - set(cols2remove_lowVars_eachPlate))

# ---------------------------------------------------
per_site_df = normalize_funcs.standardize_per_catX(
    per_site_df, "batch_plate", target_columns + uncorr_feats_cond
).copy()


control_df_perplate = (
    per_site_df.loc[per_site_df["ctrl_well"]].groupby(["batch_plate"])[target_columns].mean()
)
plates_with_controls = list(
    set(per_site_df["batch_plate"].unique().tolist())
    & set(control_df_perplate.index.unique().tolist())
)

if 1:
    per_site_df = per_site_df[per_site_df["batch_plate"].isin(plates_with_controls)].reset_index(
        drop=True
    )

    df_rep_level_scaled_meanSub = per_site_df.groupby("batch_plate")[target_columns].apply(
        subtract_control
    )

    peak_slope = np.apply_along_axis(find_end_slope2, 1, df_rep_level_scaled_meanSub.values)

#     slope = df_rep_level_scaled_meanSub.apply(
#         lambda x: find_end_slope2(x)[1], axis=1
#     )
else:
    slope = per_site_df[target_columns].apply(lambda x: find_end_slope2(x)[1], axis=1)

per_site_df[["last_peak_loc", "slope"]] = peak_slope


per_site_df = normalize_funcs.standardize_per_catX(
    per_site_df, "batch_plate", target_columns + uncorr_feats_cond + ["last_peak_loc", "slope"]
).copy()

# %% hidden=true
dataset

# %% hidden=true
# fig, ax = plt.subplots()
# im=ax.imshow(per_site_df[per_site_df['Metadata_pert_name']=='TSC2_WT'].groupby(['batch_plate','Metadata_Well'])[target_columns].mean().values)

# %% hidden=true
# per_site_df[per_site_df["Gene"].isin(hit_list)]
# per_site_df['Metadata_pert_name'].unique()

# %% hidden=true
# batches_df

# %% [markdown]
# ## 4. find PER PLATE diff and T2test between the target pattern for each pert versus controls
# - And the same for orth features
# - We calculate test stats per plate and average values across plates

# %%
dataset = "jump_compound"

# %%
# import pingouin
root_res_dir = mito_project_root_dir + "workspace/"
write_res_path = root_res_dir + "/results/virtual_screen/"
f_substr = "MeanFrac"
##############################################
from scipy.stats import norm, ttest_ind


def cohens_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(
        ((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof
    )
    return (np.mean(x) - np.mean(y)) / pooled_std


# Function to convert t-statistic to z-score
def t_to_z(t_stat, df):
    return t_stat / np.sqrt(df / (df + t_stat**2))


# Function to calculate standardized p-value from z-score
def z_to_p(z):
    return 2 * (1 - norm.cdf(abs(z)))


# sort_by_feature="Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16"
target_columns = [
    "Cells_RadialDistribution_" + f_substr + "_mito_tubeness_" + str(i) + "of16"
    for i in range(5, 17)
]

pert_col = ds_info_dict[dataset]["pert_col"]
meta_cols = ds_info_dict[dataset]["meta_cols"]


results = annot[meta_cols].drop_duplicates().reset_index(drop=True)

feature_list2 = per_site_df.columns[
    per_site_df.columns.str.contains("Cells_|Nuclei_|Cytoplasm_")
].tolist()

# list_of_cols_2add=['p_'+f for f in feature_list2]+['t_'+f for f in feature_list2]+['Count_Cells_avg']
# results=results.reindex(columns=results.columns.tolist() + list_of_cols_2add)

if "Metadata_pert_type" in per_site_df.columns:
    perts = per_site_df[per_site_df["Metadata_pert_type"].isin(["trt", "Treated"])][
        pert_col
    ].unique()

else:
    perts = per_site_df[~per_site_df["ctrl_well"]][pert_col].unique()

for peri, pert in enumerate(perts):
    if peri % 100 == 0:
        print(peri, "/", len(perts))

    per_site_df_pert = per_site_df[per_site_df[pert_col] == pert].reset_index(drop=True)
    if not per_site_df_pert.empty:
        plates_pert = (
            per_site_df_pert.groupby(["batch_plate"])
            .filter(lambda x: len(x) > 1)["batch_plate"]
            .unique()
        )
        #         plates_pert=per_site_df_pert['batch_plate'].unique()

        if len(plates_pert) > 0:
            pert_cell_count_perSite_all_plates = []
            pert_pvals_all_plates = np.full((len(plates_pert), 6), np.nan)
            pert_tvals_all_plates = np.full((len(plates_pert), 4), np.nan)
            peak_slope_all_plates = np.full((len(plates_pert), 2), np.nan)
            for pi, plate in enumerate(plates_pert):
                per_site_df_pert_plate = per_site_df_pert[
                    per_site_df_pert["batch_plate"] == plate
                ].reset_index(drop=True)
                #             if per_site_df_pert_plate.shape[0]>1:
                pert_cell_count_perSite_all_plates.append(
                    per_site_df_pert_plate["Count_Cells"].mean()
                )

                #         control_df=per_site_df[(per_site_df['pert_type']=='control') & \
                #         control_df=per_site_df[(per_site_df['Symbol'].isin(['LacZ','BFP','HcRed','LUCIFERASE'])) & \
                control_df = per_site_df[
                    (per_site_df["ctrl_well"]) & (per_site_df["batch_plate"] == plate)
                ].reset_index(drop=True)

                test_res = ttest_ind(
                    per_site_df_pert_plate["slope"], control_df["slope"], equal_var=False
                )
                cohend = cohens_d(per_site_df_pert_plate["slope"], control_df["slope"])

                pert_tvals_all_plates[pi, 3] = cohend

                degfree = (
                    per_site_df_pert_plate["slope"].shape[0] + control_df["slope"].shape[0] - 2
                )
                z_score = t_to_z(test_res.statistic, degfree)
                std_p_val = z_to_p(z_score)
                pert_pvals_all_plates[pi, 3] = std_p_val

                #         sfsdfsdfs
                #             print(test_res)
                pert_pvals_all_plates[pi, 2] = test_res.pvalue
                pert_tvals_all_plates[pi, 2] = test_res.statistic

                #         try:
                statistic, p_value, p_value_std_pattern = TwoSampleT2Test(
                    control_df[target_columns], per_site_df_pert_plate[target_columns]
                )

                pert_pvals_all_plates[pi, 0] = p_value
                pert_tvals_all_plates[pi, 0] = statistic
                pert_pvals_all_plates[pi, 4] = p_value_std_pattern

                statistic1, p_value1 = HotellingsT_internal(
                    control_df[target_columns], per_site_df_pert_plate[target_columns]
                )

                #         try:
                statistic, p_value, p_value_std_orth = TwoSampleT2Test(
                    control_df[uncorr_feats_cond], per_site_df_pert_plate[uncorr_feats_cond]
                )
                pert_pvals_all_plates[pi, 1] = p_value
                pert_tvals_all_plates[pi, 1] = statistic
                pert_pvals_all_plates[pi, 5] = p_value_std_orth

                #         except Exception as e:
                #             print("error 2")
                #             continue

                peak_slope_all_plates[pi, :] = per_site_df_pert_plate[
                    ["last_peak_loc", "slope"]
                ].median()
            #         diff_pattern=per_site_df_pert_plate[target_columns].mean()-\
            #                      control_df[target_columns].mean()
            #         peak_slope_all_plates[pi,:]=find_end_slope(diff_pattern.values)

            med_t = np.nanpercentile(pert_tvals_all_plates, 50, axis=0, interpolation="nearest")

            #     median_t_indx=[np.argwhere(pert_tvals_all_plates[:,i]==med_t[i])[0][0] if ~np.isnan(med_t[i]) else np.nan\
            #                for i in range(3)] # commneted this because I think we should stick to one metric and select the
            #                     plate based on that and then take pattern and orth t or p stats for specific plate as well
            #     if the selection is based on pattern metric then median_selection_ind=0
            median_selection_ind = 3
            if ~np.isnan(med_t[median_selection_ind]):
                median_t_indx_val = np.argwhere(
                    pert_tvals_all_plates[:, median_selection_ind] == med_t[median_selection_ind]
                )[0][0]
                median_t_indx = [median_t_indx_val] * 4
                median_t_indx_p = [median_t_indx_val] * 6
            else:
                median_t_indx = [np.nan] * 4
                median_t_indx_p = [np.nan] * 6

            results.loc[results[pert_col] == pert, "Count_Cells_avg"] = np.mean(
                pert_cell_count_perSite_all_plates
            )
            #     [pert_cell_count_perSite_all_plates[median_t_indx[i],i] for i in range(len(feature_list2))]
            results.loc[
                results[pert_col] == pert,
                [
                    "p_target_pattern",
                    "p_orth",
                    "p_slope",
                    "p_slope_std",
                    "p_pattern_std",
                    "p_orth_std",
                ],
            ] = [pert_pvals_all_plates[median_t_indx_p[i], i] for i in range(6)]
            results.loc[
                results[pert_col] == pert, ["t_target_pattern", "t_orth", "t_slope", "d_slope"]
            ] = [
                pert_tvals_all_plates[median_t_indx[i], i] if ~np.isnan(med_t[i]) else np.nan
                for i in range(4)
            ]

            results.loc[results[pert_col] == pert, ["last_peak_ind", "slope"]] = np.nanmedian(
                peak_slope_all_plates, axis=0
            )

results.sort_values(by=["slope"], ascending=False).to_csv(
    write_res_path + "/" + dataset + "_results_pattern_aug_070624.csv", index=False
)

# %%
peri

# %%
# ls /home/ubuntu/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace//results/reverse_phenotype_strength/
