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

# %%
import urllib.request

import blitzgsea as blitz
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# import mygene

# %%
# import tsquared
# import hotelling

# %%
# 'STAG1', 'ASH1L', 'ZMYM2', 'KDM6B', 'SRRM2', 'HIST1H1E'

# %%
# home_path="/home/ubuntu/" # ec2
home_path = "/home/jupyter-mhaghigh@broadinst-ee45a/"  # dgx
mito_project_root_dir = (
    home_path + "bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/"
)

# %%
import numpy as np


# import scienceplots
# plt.style.use(['science','no-latex'])
def top_table_custom(signature, library, result, set_label="Set", n=10):
    alpha_rank_lines = 0.6  # default 0.3
    sig = signature.sort_values(1, ascending=False).set_index(0)
    sig = sig[~sig.index.duplicated(keep="first")]

    plt.ioff()
    fig = plt.figure(figsize=(5, 0.5 * n), frameon=False)
    ax = fig.add_subplot(111)
    fig.patch.set_visible(False)
    plt.axis("off")

    #     ax.vlines(x=[-0.2,0,0.2,0.8], ymin=-0.1, ymax=1, color="black")
    ax.vlines(x=[0, 0.2, 0.8], ymin=-0.1, ymax=1, color="black")
    ln = np.linspace(-0.1, 1, n + 1)[::-1]
    ax.hlines(y=ln, xmin=-0.2, xmax=1, color="black")

    ax.text(-0.2 + 0.03, 1.03, "FDR", fontsize=12)
    ax.text(0.03, 1.03, "NES", fontsize=12)
    ax.text(0.84, 1.03, set_label, fontsize=12)

    for i in range(n):
        ax.text(
            -0.2 + 0.03,
            (ln[i] + ln[i + 1]) / 2,
            f"{result.iloc[i, 4]:.3f}",
            verticalalignment="center",
        )
        ax.text(
            0.03, (ln[i] + ln[i + 1]) / 2, f"{result.iloc[i, 1]:.3f}", verticalalignment="center"
        )
        ax.text(0.84, (ln[i] + ln[i + 1]) / 2, result.index[i], verticalalignment="center")

        gs = set(library[result.index[i]])
        hits = np.array([i for i, x in enumerate(sig.index) if x in gs])
        hits = (hits / len(sig.index)) * 0.6 + 0.2

        if result.iloc[i, 1] > 0:
            ax.vlines(hits, ymax=ln[i], ymin=ln[i + 1], color="red", lw=0.5, alpha=alpha_rank_lines)
        else:
            ax.vlines(
                hits, ymax=ln[i], ymin=ln[i + 1], color="blue", lw=0.5, alpha=alpha_rank_lines
            )
    fig.patch.set_facecolor("white")
    plt.ion()
    return fig


# p_val_orth_col="p_orth_std";#p_orth_std
orth_bh_corrected_critical_dict = {
    "taorf": 0.0328682154922577,
    "lincs": 0.0482929972989429,
    "CDRP": 0.0481242427803521,
    "jump_orf": 0.0375156321477045,
    "jump_crispr": 0.0478530069113243,
    "jump_compound": 0.0474630677339327,
}
datasets_info_dict = {
    "taorf": {"key_col": "Metadata_gene_name", "key_sample_id_col": "Metadata_broad_sample"},
    "jump_orf": {"key_col": "Metadata_Symbol", "key_sample_id_col": "Metadata_broad_sample"},
    "jump_crispr": {"key_col": "Metadata_Symbol", "key_sample_id_col": "Metadata_JCP2022"},
    "jump_compound": {"key_col": "Metadata_Symbol", "key_col_ref_set": "Metadata_JCP2022"},
    "lincs": {"key_col": "Metadata_pert_name", "key_col_ref_set": "Metadata_pert_id_dose"},
    "CDRP": {"key_col": "Metadata_pert_id", "key_col_ref_set": "Metadata_Sample_Dose"},
}

# Metadata_Sample_Dose
# key_col_ref_set='Metadata_Sample_Dose'
# Metadata_pert_name_lowercase, Metadata_pert_id_dose

# %%

# %%

# %% [markdown] heading_collapsed=true
# ## Set Enrichment Analysis was done using blitzGSEA package
# - https://github.com/MaayanLab/blitzgsea

# %% hidden=true
# # ls ~/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/results/virtual_screen

# %% [markdown] heading_collapsed=true
# ## GO - ORF datasets

# %% hidden=true
fName = "d_slope"

write_res_path = mito_project_root_dir + "workspace/results/virtual_screen/"


cell_count_filter_enabled = True
orth_filter_enabled = False


database_ls = [
    "OMIM_Expanded",
    "OMIM_Disease",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "Human_Phenotype_Ontology",
    "MGI_Mammalian_Phenotype_Level_4_2021",
    "KEGG_2021_Human",
    "WikiPathway_2021_Human",
    "PFOCR_Pathways_2023",
    "Proteomics_Drug_Atlas_2023",
    "GWAS_Catalog_2023",
]

plt.close("all")
for dataset in ["jump_crispr"]:  # ['taorf','jump_orf','jump_crispr']:
    key_col = datasets_info_dict[dataset]["key_col"]
    res_df_jumporf0 = pd.read_csv(write_res_path + dataset + "_results_pattern_aug_070624.csv")
    # res_df_jumporf0['t_target_pattern_signed']=np.sign(res_df_jumporf0['t_slope'])*res_df_jumporf0['t_target_pattern']
    res_df_jumporf0 = res_df_jumporf0[
        ~(res_df_jumporf0[key_col].isnull()) & ~(res_df_jumporf0[fName].isnull())
    ].reset_index(drop=True)
    res_df_jumporf = res_df_jumporf0.groupby(key_col).median(numeric_only=True).reset_index()

    # print(res_df_jumporf.shape)
    if cell_count_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["Count_Cells_avg"] > res_df_jumporf["Count_Cells_avg"].quantile(0.1)
        ].reset_index(drop=True)

    if orth_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["p_orth_std"] > orth_bh_corrected_critical_dict[dataset]
        ].reset_index(drop=True)
    # print(res_df_jumporf.shape)

    sig_df = (
        res_df_jumporf[[key_col, fName]]
        .rename(columns={key_col: 0, fName: 1})
        .sort_values(by=1)
        .reset_index(drop=True)
    )

    for database_str in database_ls:
        library = blitz.enrichr.get_library(database_str)
        #         sig_df[1]=sig_df[1]/10
        result_orf = blitz.gsea(sig_df, library, min_size=4, seed=1)

        top_res_df = result_orf[result_orf["fdr"] < 0.05]
        if len(top_res_df) > 0:
            # fig_table = top_table_custom(sig_df, library, top_res_df,set_label="Gene Set: "+database_str, n=top_res_df.shape[0])
            fig_table = top_table_custom(
                sig_df, library, top_res_df, set_label="Gene Set ", n=top_res_df.shape[0]
            )
            plt.show()
        #     fig_table.savefig("top_10_"+dataset+"_"+database_str+".png", bbox_inches='tight')
        if top_res_df.shape[0] > 0:
            print(
                " * ",
                dataset,
                ": Out of ",
                result_orf.shape[0],
                database_str,
                " categories, ",
                top_res_df.shape[0],
                " are enriched.",
            )
            # print(top_res_df[['nes','fdr','leading_edge']].to_markdown())
#             fig_table = top_table_custom(sig_df, library, result_orf,set_label="Gene Set", n=top_res_df.shape[0])
# fig_table.savefig(dataset+"_top_table.png", bbox_inches='tight')
#             plt.show()
#         result_orf.head(10)


# %%
# top_res_df

# %% hidden=true
# # res_df_jumporf0 = pd.read_csv(write_res_path+dataset+'_results_pattern_1.csv')
# res_df_jumporf0=pd.read_csv(write_res_path+dataset+'_results_pattern_aug_070624.csv')
# res_df_jumporf0.groupby(key_col).median().reset_index()

# %% [markdown] heading_collapsed=true
# ## Table 3 Genes

# %% hidden=true
fName = "d_slope"
write_res_path = mito_project_root_dir + "workspace/results/virtual_screen/"


cell_count_filter_enabled = False
orth_filter_enabled = False

# ------------------------------

# Fusion Fission Transport Mitophagy
Fusion_ls = [
    "AFG3L2",
    "BAK1",
    "BAX",
    "BNIP3",
    "CHCHD3",
    "DNM1L",
    "FIS1",
    "GDAP1",
    "MFF",
    "MFN1",
    "MFN2",
    "MIEF1",
    "MIEF2",
    "MUL1",
    "OMA1",
    "OPA1",
    "PLD6",
    "PRKN",
    "STOML2",
    "USP30",
]
Fission_ls = [
    "BNIP3",
    "COX10",
    "DHODH",
    "DNM1L",
    "DNM2",
    "DNM3",
    "FIS1",
    "GDAP1",
    "LRRK2",
    "MARCH5",
    "MFF",
    "MIEF1",
    "MIEF2",
    "MTFP1",
    "MTFR1",
    "MTFR1L",
    "MTFR2",
    "MUL1",
    "MYO19",
    "OPA1",
    "PINK1",
    "PRKN",
    "SLC25A46",
]
Transport_ls = [
    "CLUH",
    "DNM1L",
    "KIF1B",
    "LRPPRC",
    "LRRK2",
    "MAIP1",
    "MAP1S",
    "MFN1",
    "MFN2",
    "MGARP",
    "MST01",
    "MUL1",
    "OPA1",
    "RHOT1",
    "RHOT2",
    "SYNJ2BP",
    "TRAK1",
    "TRAK2",
    "UBB",
]
Mitophagy_ls = [
    "ATPIF1",
    "BNIP3",
    "BNIP3L",
    "CISD2",
    "DNM1L",
    "FIS1",
    "FUNDC1",
    "FUNDC2",
    "HK2",
    "HTRA2",
    "MFN2",
    "MUL1",
    "PARK7",
    "PINK1",
    "SQSTM1",
    "TOMM7",
    "TSPO",
    "VDAC1",
]


table3_list = [
    "AFG3L2",
    "BAK1",
    "BAX",
    "BNIP3",
    "CHCHD3",
    "DNM1L",
    "FIS1",
    "GDAP1",
    "MFF",
    "MFN1",
    "MFN2",
    "MIEF1",
    "MIEF2",
    "MUL1",
    "OMA1",
    "OPA1",
    "PLD6",
    "PRKN",
    "STOML2",
    "USP30",
    "BNIP3",
    "COX10",
    "DHODH",
    "DNM1L",
    "DNM2",
    "DNM3",
    "FIS1",
    "GDAP1",
    "LRRK2",
    "MARCH5",
    "MFF",
    "MIEF1",
    "MIEF2",
    "MTFP1",
    "MTFR1",
    "MTFR1L",
    "MTFR2",
    "MUL1",
    "MYO19",
    "OPA1",
    "PINK1",
    "PRKN",
    "SLC25A46",
    "CLUH",
    "DNM1L",
    "KIF1B",
    "LRPPRC",
    "LRRK2",
    "MAIP1",
    "MAP1S",
    "MFN1",
    "MFN2",
    "MGARP",
    "MST01",
    "MUL1",
    "OPA1",
    "RHOT1",
    "RHOT2",
    "SYNJ2BP",
    "TRAK1",
    "TRAK2",
    "UBB",
    "ATPIF1",
    "BNIP3",
    "BNIP3L",
    "CISD2",
    "DNM1L",
    "FIS1",
    "FUNDC1",
    "FUNDC2",
    "HK2",
    "HTRA2",
    "MFN2",
    "MUL1",
    "PARK7",
    "PINK1",
    "SQSTM1",
    "TOMM7",
    "TSPO",
    "VDAC1",
]
# library_table3={'Table3genes':list(set(table3_list)),'Fusion':Fusion_ls,'Fission':Fission_ls,\
#                 'Transport':Transport_ls,'Mitophagy':Mitophagy_ls}

# key_sample_id_col='Metadata_broad_sample'
plt.close("all")
for dataset in ["jump_orf", "jump_crispr"]:
    key_col = datasets_info_dict[dataset]["key_col"]
    #     key_sample_id_col=datasets_info_dict[dataset]['key_sample_id_col'];
    key_sample_id_col = datasets_info_dict[dataset]["key_sample_id_col"]
    res_df_jumporf0 = pd.read_csv(write_res_path + dataset + "_results_pattern_aug_070624.csv")
    res_df_jumporf0["t_target_pattern_signed"] = (
        np.sign(res_df_jumporf0["t_slope"]) * res_df_jumporf0["t_target_pattern"]
    )
    res_df_jumporf0 = res_df_jumporf0[
        ~(res_df_jumporf0[key_col].isnull()) & ~(res_df_jumporf0[fName].isnull())
    ].reset_index(drop=True)
    #     res_df_jumporf = res_df_jumporf0.groupby(key_col).median(numeric_only=True).reset_index()

    # Define the aggregation dictionary
    agg_dict = {
        col: "median" if pd.api.types.is_numeric_dtype(res_df_jumporf0[col]) else "first"
        for col in res_df_jumporf0.columns
        if col != key_sample_id_col
    }

    # Group by key_sample_id_col and aggregate
    res_df_jumporf = res_df_jumporf0.groupby(key_sample_id_col).agg(agg_dict).reset_index()

    Transport_ls_brd = (
        res_df_jumporf.loc[res_df_jumporf[key_col].isin(Transport_ls), key_sample_id_col]
        .unique()
        .tolist()
    )

    table3_ls_brd = (
        res_df_jumporf.loc[res_df_jumporf[key_col].isin(table3_list), key_sample_id_col]
        .unique()
        .tolist()
    )

    # library_table3={'Table3genes':list(set(table3_list)),'Transport':Transport_ls}
    # library_table3={'Table3genes':list(set(table3_list)),'Transport':Transport_ls}
    library_table3 = {
        "Mitochondrial dynamics": list(set(table3_ls_brd)),
        "Mitochondrial transport": Transport_ls_brd,
    }

    # print(res_df_jumporf.shape)
    if cell_count_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["Count_Cells_avg"] > res_df_jumporf["Count_Cells_avg"].quantile(0.1)
        ].reset_index(drop=True)

    if orth_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["p_orth_std"] > orth_bh_corrected_critical_dict[dataset]
        ].reset_index(drop=True)
    # print(res_df_jumporf.shape)

    sig_df = (
        res_df_jumporf[[key_sample_id_col, fName]]
        .rename(columns={key_sample_id_col: 0, fName: 1})
        .sort_values(by=1)
        .reset_index(drop=True)
    )

    result_table3set = blitz.gsea(sig_df, library_table3, seed=1)

    print(
        result_table3set[["es", "nes", "pval", "sidak", "fdr", "geneset_size"]].to_markdown(
            index=True
        )
    )

    fig_table = top_table_custom(
        sig_df, library_table3, result_table3set.sort_values(by="Term"), set_label="Gene Set", n=2
    )
    # fig_table.savefig(dataset+"_top_table.png", bbox_inches='tight')
    plt.show()
    result_table3set


# %%

# %% hidden=true
# agg_dict

# %% hidden=true
print(
    result_table3set[["es", "nes", "pval", "sidak", "fdr", "geneset_size"]].to_markdown(index=True)
)

# %% hidden=true

# %% hidden=true
mito_project_root_dir

# %% [markdown]
# ## SZ related genes according to "Rare coding variants in ten genes confer substantial risk for schizophrenia" nature paper

# %%
fName = "d_slope"
write_res_path = mito_project_root_dir + "workspace/results/virtual_screen/"


cell_count_filter_enabled = True
orth_filter_enabled = False

# ------------------------------

SZ_genes_Daly = [
    "CACNA1G",
    "GRIN2A",
    "GRIA3",
    "TRIO",
    "SP4",
    "RB1CC1",
    "SETD1A",
    "XPO7",
    "CUL1",
    "HERC1",
]
DD_genes_Daly = ["STAG1", "ASH1L", "ZMYM2", "KDM6B", "SRRM2", "HIST1H1E", "PMEPA1"]
misc_genes = ["CACNA1A", "AKAP11", "SETD1A"]


# key_sample_id_col='Metadata_broad_sample'

plt.close("all")
for dataset in ["jump_orf", "jump_crispr"]:
    key_col = datasets_info_dict[dataset]["key_col"]
    #     key_sample_id_col=datasets_info_dict[dataset]['key_sample_id_col'];
    key_sample_id_col = datasets_info_dict[dataset]["key_sample_id_col"]
    res_df_jumporf0 = pd.read_csv(write_res_path + dataset + "_results_pattern_aug_070624.csv")
    res_df_jumporf0 = res_df_jumporf0[
        ~(res_df_jumporf0[key_col].isnull()) & ~(res_df_jumporf0[fName].isnull())
    ].reset_index(drop=True)
    #     res_df_jumporf = res_df_jumporf0.groupby(key_col).median(numeric_only=True).reset_index()

    # Define the aggregation dictionary
    agg_dict = {
        col: "median" if pd.api.types.is_numeric_dtype(res_df_jumporf0[col]) else "first"
        for col in res_df_jumporf0.columns
        if col != key_sample_id_col
    }

    # Group by key_sample_id_col and aggregate
    res_df_jumporf = res_df_jumporf0.groupby(key_sample_id_col).agg(agg_dict).reset_index()

    library_Daly = {}
    library_Daly_symbols = {}
    gls_names = ["SZ genes", "DD/ID genes", "Misc genes"]
    for gls_name, gls in zip(gls_names, [SZ_genes_Daly, DD_genes_Daly, misc_genes], strict=False):
        library_Daly[gls_name] = list(
            set(
                res_df_jumporf.loc[res_df_jumporf[key_col].isin(gls), key_sample_id_col]
                .unique()
                .tolist()
            )
        )
        library_Daly_symbols[gls_name] = list(
            set(res_df_jumporf.loc[res_df_jumporf[key_col].isin(gls), key_col].unique().tolist())
        )

    print("## ", dataset)
    print(library_Daly_symbols)

    # print(res_df_jumporf.shape)
    if cell_count_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["Count_Cells_avg"] > res_df_jumporf["Count_Cells_avg"].quantile(0.1)
        ].reset_index(drop=True)

    if orth_filter_enabled:
        res_df_jumporf = res_df_jumporf[
            res_df_jumporf["p_orth_std"] > orth_bh_corrected_critical_dict[dataset]
        ].reset_index(drop=True)
    # print(res_df_jumporf.shape)

    sig_df = (
        res_df_jumporf[[key_sample_id_col, fName]]
        .rename(columns={key_sample_id_col: 0, fName: 1})
        .sort_values(by=1)
        .reset_index(drop=True)
    )

    result_Daly = blitz.gsea(sig_df, library_Daly, min_size=1, seed=1)

    print(
        result_Daly[["es", "nes", "pval", "sidak", "fdr", "geneset_size"]].to_markdown(index=True)
    )

    fig_table = top_table_custom(sig_df, library_Daly, result_Daly, set_label="Gene Set", n=2)
    # # fig_table.savefig(dataset+"_top_table.png", bbox_inches='tight')
    plt.show()
    # result_Daly


# %%

# %% [markdown] heading_collapsed=true
# ## Rakesh Categories - LINCS

# %% [markdown] heading_collapsed=true
# ## Experiment our list towards Rakesh categories on LINCS set

# %% hidden=true
### set params

target_feat = "d_slope"


dataset = "lincs"

cell_count_filter_enabled = False
orth_filter_enabled = False

## compound selection strategies
# 1. keeping all doses


# 2. keep just the dose with maximum phenotype strength


res_df_lincs = pd.read_csv(write_res_path + dataset + "_results_pattern_aug_070624.csv")
# res_df_lincs['t_target_pattern_signed']=np.sign(res_df_lincs['t_slope'])*res_df_lincs['t_target_pattern']
# res_df_lincs[target_feat+'_abs']=res_df_lincs[target_feat].abs()

res_df_lincs = res_df_lincs[
    ~(res_df_lincs[target_feat].isnull()) & ~(res_df_lincs["Metadata_pert_name"].isnull())
].reset_index(drop=True)

res_df_lincs[target_feat + "_abs"] = res_df_lincs[target_feat].abs()

drug_list_rakesh0 = pd.read_excel(
    mito_project_root_dir + "workspace/metadata/CompoundClusters08202104.xlsx"
)


# 'lincs':{'key_col':'Metadata_pert_name',\
# 'key_col_ref_set':'Metadata_pert_id_dose'},\

if cell_count_filter_enabled:
    res_df_lincs = res_df_lincs[
        res_df_lincs["Count_Cells_avg"] > res_df_lincs["Count_Cells_avg"].quantile(0.1)
    ].reset_index(drop=True)
if orth_filter_enabled:
    res_df_lincs = res_df_lincs[
        res_df_lincs["p_orth_std"] > orth_bh_corrected_critical_dict[dataset]
    ].reset_index(drop=True)


lincs_prt_names = set(res_df_lincs["Metadata_pert_name"].str.lower().unique())
res_df_lincs["Metadata_pert_name_lowercase"] = res_df_lincs["Metadata_pert_name"].str.lower()

if 1:
    abs_max_indices = (
        res_df_lincs.groupby(["Metadata_pert_id"])[target_feat + "_abs"].idxmax().values
    )
    abs_max_indices = (
        res_df_lincs.groupby(["Metadata_pert_name_lowercase"])[target_feat + "_abs"].idxmax().values
    )
    res_df_lincs = res_df_lincs.loc[abs_max_indices].reset_index(drop=True)


rakesh_cats = drug_list_rakesh0.columns.tolist()
# ap_table=pd.DataFrame(index=rakesh_cats,columns=['AP','AP-abs','Baseline-AP'])

key_col_ref_set = "Metadata_pert_name_lowercase"
key_col_ref_set = "Metadata_pert_id"
# key_col_ref_set='Metadata_pert_id_dose'


lib_rakesh = {}
for r in rakesh_cats:
    overlap_drugs = set(drug_list_rakesh0[r].str.lower().unique()) & lincs_prt_names
    overlap_drugs = {x for x in overlap_drugs if x == x}
    if len(overlap_drugs) > 0:
        #         lib_rakesh[r]=list(overlap_drugs)
        ## if adding all doses
        #         pert_id_dose_ls=res_df_lincs2.loc[res_df_lincs2['Metadata_pert_name_lowercase'].isin(overlap_drugs),'Metadata_pert_name_lowercase'].tolist()
        pert_id_dose_ls = (
            res_df_lincs.loc[
                res_df_lincs["Metadata_pert_name_lowercase"].isin(overlap_drugs), key_col_ref_set
            ]
            .unique()
            .tolist()
        )
        if len(pert_id_dose_ls) > 1:
            lib_rakesh[r] = list(pert_id_dose_ls)

    print(
        "-",
        r,
        len(overlap_drugs),
        " out of ",
        len(set(drug_list_rakesh0[r].str.lower().unique())),
        " exists in LINCS set.",
    )
####################################


sig_df_lincs = (
    res_df_lincs[[key_col_ref_set, target_feat]]
    .rename(columns={key_col_ref_set: 0, target_feat: 1})
    .sort_values(by=1)
    .reset_index(drop=True)
)


# sig_df_lincs=sig_df_lincs.groupby([0]).mean().reset_index()
result_rakesh_enrich = blitz.gsea(sig_df_lincs, lib_rakesh, min_size=5, seed=1)
plt.close("all")
####################################
fig_table = top_table_custom(
    sig_df_lincs,
    lib_rakesh,
    result_rakesh_enrich,
    set_label="Drug Set",
    n=result_rakesh_enrich.shape[0],
)
plt.show()
# fig_table.savefig("ES_rakesh_new.png", bbox_inches='tight')

# %%
result_rakesh_enrich

# %%
lib_rakesh.keys()

# %% hidden=true
# sig_df_lincs.groupby([0]).mean().reset_index()
# res_df_lincs.groupby('Metadata_pert_name_lowercase').size().sort_values()

# %% hidden=true
# sig_df_lincs.groupby([0]).size().sort_values()
# sig_df_lincs

# %% hidden=true
print(
    result_rakesh_enrich.rename(columns={"geneset_size": "drugset_size"})[
        ["es", "nes", "pval", "fdr", "drugset_size"]
    ].to_markdown(index=True)
)

# %% hidden=true
result_rakesh_enrich.columns

# %% hidden=true
fig_table = top_table_custom(sig_df_lincs, lib_rakesh, result_rakesh_enrich, n=8)
fig_table.savefig("ES_rakesh.png", bbox_inches="tight")

# %% hidden=true
# result_rakesh_enrich = blitz.gsea(sig_df_lincs, lib_rakesh)

# # fig = blitz.plot.running_sum(sig_df_lincs, "Table3genes", lib_rakesh, result=result_rakesh_enrich, compact=False)
# # fig.savefig("ES_rakesh1.png", bbox_inches='tight')

# # fig_compact = blitz.plot.running_sum(sig_df_lincs, "Table3genes", lib_rakesh, result=result_rakesh_enrich, compact=True)
# # fig_compact.savefig("ES_rakesh2.png", bbox_inches='tight')

# fig_table = blitz.plot.top_table(sig_df_lincs, lib_rakesh, result_rakesh_enrich, n=5)
# fig_table.savefig("ES_rakesh_max.png", bbox_inches='tight')

fig_table = top_table_custom(sig_df_lincs, lib_rakesh, result_rakesh_enrich, n=5, seed=1)
fig_table.savefig("ES_rakesh_max2.png", bbox_inches="tight")


# %% hidden=true
for r in lib_rakesh:
    fig = blitz.plot.running_sum(
        sig_df_lincs, r, lib_rakesh, result=result_rakesh_enrich, compact=False
    )
    fig.savefig("ES_rakesh_" + r + ".png", bbox_inches="tight")

# %% hidden=true
print(result_rakesh_enrich.rename(columns={"geneset_size": "drugset_size"}).to_markdown(index=True))

# %% hidden=true
result_rakesh_enrich

# %% hidden=true
# res_df_lincs[res_df_lincs['Metadata_pert_name_lowercase']=='sertraline']

# %% hidden=true
# res_df_lincs['Metadata_moa'].shape,
# res_df_lincs['Metadata_moa'].unique()

# %% hidden=true

# %% hidden=true

# %% hidden=true
# res_df_lincs_nonnan_moa=res_df_lincs[~res_df_lincs['Metadata_moa'].isnull()].reset_index(drop=True)

# %% hidden=true
# res_df_lincs_nonnan_moa.loc[res_df_lincs_nonnan_moa['Metadata_moa'].str.contains('norepinephrine reuptake inhibitor'),'Metadata_moa'].unique()

# %% [markdown]
# ## Enrichement of MOA classes - new pattern

# %%
target_feat = "d_slope"
# target_feat="t_target_pattern_signed"

write_res_path = mito_project_root_dir + "workspace/results/virtual_screen/"

cell_count_filter_enabled = True
orth_filter_enabled = False


for dataset in ["lincs", "CDRP"]:  # ['lincs','CDRP']:
    key_col = datasets_info_dict[dataset]["key_col"]
    key_col_ref_set = datasets_info_dict[dataset]["key_col_ref_set"]

    # res_df_lincs = pd.read_csv(write_res_path+dataset+'_results_pattern_2.csv')
    #     res_df_lincs = pd.read_csv(write_res_path+dataset+'_results_pattern_aug.csv')
    res_df_lincs = pd.read_csv(write_res_path + dataset + "_results_pattern_aug_070624.csv")
    res_df_lincs["t_target_pattern_signed"] = (
        np.sign(res_df_lincs["t_slope"]) * res_df_lincs["t_target_pattern"]
    )

    res_df_lincs["Metadata_moa"] = res_df_lincs["Metadata_moa"].str.lower()

    if cell_count_filter_enabled:
        res_df_lincs = res_df_lincs[
            res_df_lincs["Count_Cells_avg"] > res_df_lincs["Count_Cells_avg"].quantile(0.1)
        ].reset_index(drop=True)

    if orth_filter_enabled:
        res_df_lincs = res_df_lincs[
            res_df_lincs["p_orth_std"] > orth_bh_corrected_critical_dict[dataset]
        ].reset_index(drop=True)

    res_df_lincs[target_feat + "_abs"] = res_df_lincs[target_feat].abs()
    res_df_lincs["t_orth_abs"] = res_df_lincs["t_orth"].abs()

    res_df_lincs = res_df_lincs[
        ~(res_df_lincs[target_feat].isnull()) & ~(res_df_lincs[key_col].isnull())
    ].reset_index(drop=True)
    res_df_lincs[key_col + "_lowercase"] = res_df_lincs[key_col].str.lower()

    if 0:
        abs_max_indices = (
            res_df_lincs.groupby([key_col + "_lowercase"])[target_feat + "_abs"].idxmax().values
        )
        #         abs_max_indices=res_df_lincs.groupby([key_col+'_lowercase'])["t_orth_abs"].idxmin().values
        # abs_max_indices=res_df_lincs.groupby([key_col+'_lowercase'])["p_orth_std"].idxmax().values
        res_df_lincs = res_df_lincs.loc[abs_max_indices].reset_index(drop=True)

    if 0:
        res_df_lincs_nonnan_moa = res_df_lincs[~res_df_lincs["Metadata_moa"].isnull()].reset_index(
            drop=True
        )
    else:
        res_df_lincs_nonnan_moa = res_df_lincs.copy()
        res_df_lincs_nonnan_moa["Metadata_moa"] = res_df_lincs_nonnan_moa["Metadata_moa"].fillna(
            "nan"
        )

    # res_df_lincs_nonnan_moa = res_df_lincs[~res_df_lincs['Metadata_moa'].isnull()].reset_index(drop=True)

    # key_col_ref_set='Metadata_pert_name_lowercase'
    #     key_col_ref_set='Metadata_pert_id_dose'

    sig_df_lincs_2 = (
        res_df_lincs_nonnan_moa[[key_col_ref_set, target_feat]]
        .rename(columns={key_col_ref_set: 0, target_feat: 1})
        .sort_values(by=1)
        .reset_index(drop=True)
    )

    unq_moas0 = []
    all_moas = res_df_lincs["Metadata_moa"].unique().tolist()

    for mi in all_moas:
        if mi == mi:
            unq_moas0 = unq_moas0 + mi.split("|")

    unq_moas = list(set(unq_moas0))  # [50:60]

    lib_moa = {}
    for r in unq_moas:
        pert_id_dose_ls = (
            res_df_lincs_nonnan_moa.loc[
                res_df_lincs_nonnan_moa["Metadata_moa"].str.contains(r), key_col_ref_set
            ]
            .unique()
            .tolist()
        )
        if len(pert_id_dose_ls) > 1:
            lib_moa[r] = list(pert_id_dose_ls)

    min_n_compounds_per_cat = 10
    result_moa_enrich = blitz.gsea(
        sig_df_lincs_2, lib_moa, min_size=min_n_compounds_per_cat, seed=1
    )
    enriched_df = result_moa_enrich[result_moa_enrich["fdr"] < 0.05]

    if 1:
        fig_table = top_table_custom(
            sig_df_lincs_2, lib_moa, result_moa_enrich, set_label="MOA", n=enriched_df.shape[0]
        )
        plt.show()

    print(
        dataset,
        " has ",
        str(len(unq_moas)),
        " MOA categories! ",
        str(result_moa_enrich.shape[0]),
        " have min size more than ",
        min_n_compounds_per_cat,
    )
    print(enriched_df.shape[0], "/", str(result_moa_enrich.shape[0]), " categories are enriched!")
    print(enriched_df.index)


# %%
# print(result_moa_enrich.rename(columns={'geneset_size':"moa_set_size"})[['es','nes','pval','fdr','moa_set_size']].to_markdown(index=True))

# %%
# res_df_lincs['d_slope']

# %%
# res_df_lincs_nonnan_moa.groupby('Metadata_moa').mean().sort_values(by='slope_abs').loc['mTOR inhibitor']

# %%
# fig_table = top_table_custom(sig_df_lincs_2, lib_moa, result_moa_enrich,set_label="MOA", n=20)
# plt.show()
# fig_table.savefig("ES_moa_slope_top10.png", bbox_inches='tight')

# %%
# result_moa_enrich = blitz.gsea(sig_df_lincs_2, lib_moa,min_size=2)

# %% hidden=true

# %% hidden=true
