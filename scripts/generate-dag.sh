#!/usr/bin/env bash
# Generate DAG visualizations (simplified rules + full jobs for all three methods)

set -e  # Exit on error

mkdir -p docs/pipeline

# Generate simplified DAGs (rules only)
pixi run snakemake all_baseline --rulegraph 2>&1 | tail -n +2 > /tmp/rg_baseline.dot
pixi run snakemake all_notebook --rulegraph 2>&1 | tail -n +2 > /tmp/rg_notebook.dot
pixi run snakemake all_module_csvs all_module_diagnostics --rulegraph 2>&1 | tail -n +2 > /tmp/rg_module.dot
dot -Tpng /tmp/rg_baseline.dot -o docs/pipeline/dag_baseline_rules.png
dot -Tpng /tmp/rg_notebook.dot -o docs/pipeline/dag_notebook_rules.png
dot -Tpng /tmp/rg_module.dot -o docs/pipeline/dag_module_rules.png

# Generate full DAGs (all job instances)
pixi run snakemake all_baseline --dag 2>&1 | tail -n +2 > /tmp/dag_baseline.dot
pixi run snakemake all_notebook --dag 2>&1 | tail -n +2 > /tmp/dag_notebook.dot
pixi run snakemake all_module_csvs all_module_diagnostics --dag 2>&1 | tail -n +2 > /tmp/dag_module.dot
dot -Tpng /tmp/dag_baseline.dot -o docs/pipeline/dag_baseline_jobs.png
dot -Tpng /tmp/dag_notebook.dot -o docs/pipeline/dag_notebook_jobs.png
dot -Tpng /tmp/dag_module.dot -o docs/pipeline/dag_module_jobs.png

# Cleanup temp files
rm /tmp/rg_baseline.dot /tmp/rg_notebook.dot /tmp/rg_module.dot /tmp/dag_baseline.dot /tmp/dag_notebook.dot /tmp/dag_module.dot

echo "Generated in docs/pipeline/:"
echo "  METHOD 0 (Baseline - Validated S3 CSVs):"
echo "    - dag_baseline_rules.png (simplified)"
echo "    - dag_baseline_jobs.png (full)"
echo "  METHOD 1 (Notebook - Complete Regenerated):"
echo "    - dag_notebook_rules.png (simplified)"
echo "    - dag_notebook_jobs.png (full)"
echo "  METHOD 2 (Module - Clean Regenerated, Incomplete):"
echo "    - dag_module_rules.png (simplified)"
echo "    - dag_module_jobs.png (full)"
