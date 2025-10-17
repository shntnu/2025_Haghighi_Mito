# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Test singlecell-morph import

# %%
from singlecell.preprocess import extract_cpfeature_names, handle_nans
from singlecell.process import precision_recall, statistical_tests
from singlecell.read import read_single_cell_sql
from singlecell.visualize import cluster, visualize_n_SingleCell

# %%
print("Successfully imported singlecell-morph!")
