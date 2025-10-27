#!/usr/bin/env bash
# Quick hypothesis testing loop for baseline reproduction
#
# Usage: scripts/test-hypothesis.sh [DATASET]
#   DATASET: Dataset to test (default: taorf)
#
# This script:
# 1. Removes cached CSV to force regeneration
# 2. Runs virtual screen generation
# 3. Runs diagnostics and extracts correlations
#
# Total time: ~7 seconds for taorf

set -euo pipefail

DATASET="${1:-taorf}"

echo "Testing baseline reproduction for dataset: ${DATASET}"
echo "================================================"
echo ""

# Remove cached outputs to force full regeneration
echo "Removing cached outputs..."
rm -f "data/processed/virtual_screen_module/${DATASET}_results_pattern_aug_070624.csv"
rm -f "data/processed/virtual_screen_module/${DATASET}_baseline_comparison.csv"
rm -f "data/processed/figures/diagnostics/${DATASET}_comparison_metrics.png"

# Run the tight loop
echo "Running virtual screen + diagnostics..."
echo ""

just generate-module-csv-for "${DATASET}" && \
    just diagnose-for "${DATASET}" 2>&1 | grep -E "^2025.*: r="

echo ""
echo "Done! Check correlations above."
echo "Key metrics: slope (main target), t_target_pattern (validates preprocessing)"
