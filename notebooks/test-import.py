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
from singlecell.read import read_single_cell_sql
from singlecell.preprocess import handle_nans, extract_cpfeature_names
from singlecell.visualize import visualize_n_SingleCell, cluster
from singlecell.process import statistical_tests, precision_recall

# %%
print("Successfully imported singlecell-morph!")
