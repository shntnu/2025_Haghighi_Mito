#!/usr/bin/env bash
# Generate DAG visualizations (simplified rules + full jobs for both pipelines)

set -e  # Exit on error

mkdir -p docs/pipeline

# Generate simplified DAGs (rules only)
pixi run snakemake all_baseline --rulegraph 2>&1 | tail -n +2 > /tmp/rg_baseline.dot
pixi run snakemake all_regenerated --rulegraph 2>&1 | tail -n +2 > /tmp/rg_regenerated.dot
dot -Tpng /tmp/rg_baseline.dot -o docs/pipeline/dag_baseline_rules.png
dot -Tpng /tmp/rg_regenerated.dot -o docs/pipeline/dag_regenerated_rules.png

# Generate full DAGs (all job instances)
pixi run snakemake all_baseline --dag 2>&1 | tail -n +2 > /tmp/dag_baseline.dot
pixi run snakemake all_regenerated --dag 2>&1 | tail -n +2 > /tmp/dag_regenerated.dot
dot -Tpng /tmp/dag_baseline.dot -o docs/pipeline/dag_baseline_jobs.png
dot -Tpng /tmp/dag_regenerated.dot -o docs/pipeline/dag_regenerated_jobs.png

# Cleanup temp files
rm /tmp/rg_baseline.dot /tmp/rg_regenerated.dot /tmp/dag_baseline.dot /tmp/dag_regenerated.dot

echo "Generated in docs/pipeline/:"
echo "  - dag_baseline_rules.png (simplified)"
echo "  - dag_regenerated_rules.png (simplified)"
echo "  - dag_baseline_jobs.png (full)"
echo "  - dag_regenerated_jobs.png (full)"
