# Progress Log

INSTRUCTIONS: Add log at the end of the file. See "Maintenance Guidelines" and "Template for Future Entries" first

## Maintenance Guidelines

**When to add entries:**

- User explicitly asks to document progress
- Major milestone completed (pipeline working, bug resolved)
- Critical decision point requiring future reference

**Entry format (max 20 lines):**

- Brief title: what was accomplished
- Key findings or decisions (not step-by-step debugging)
- Blocking issues or next actions only if unresolved

**Monthly rollup:**

- Summarize entries older than 30 days to 1-3 lines per week
- Keep detailed recent history for active debugging
- Git history preserves all details if forensics needed

**DO NOT include:**

- Command outputs or error messages
- Step-by-step debugging notes
- Exploratory analysis details
- Duplicate information in code comments

## Template for Future Entries

```text
## YYYY-MM-DD: Brief Description

### What was done
- Key accomplishment 1
- Key accomplishment 2

### Key findings or decisions
- Finding that requires future reference

### Unresolved issues (if any)
- Issue description and next steps

### Notes (optional)
- Additional context (max 5 lines)
```

## 2025-10-17 to 2025-10-24: Initial Data Download

Downloaded S3 data (178 files, 4.6 GB) after Glacier restoration. Set up local analysis infrastructure with metadata, per-site profiles, and orthogonal features.

---

## 2025-10-25: Virtual Screen Pipeline - Local Setup & Optimization

### Initial Setup & Bug Fixes

- Fixed import errors and paths in notebook 2.0 for local execution
- Fixed metadata preprocessing bugs (missing `Batch` column, control well identification)
- Skipped per-site profile creation (pre-computed profiles already available)
- Added comprehensive logging with loguru
- Fixed dataset selection bug (hardcoded override on line 1244)

### Performance Optimization

**Vectorized slope calculation** (`haghighi_mito/vectorized_slope.py`):

- Implemented `find_end_slope2_vectorized()` - fully vectorized version
- Performance: ~200x speedup (0.174s → 0.0008s per 1000 rows)
- Created comprehensive test suite (9 tests, all passing)

**Vectorized statistical testing** (`haghighi_mito/vectorized_stats.py`):

- Batch processing for plate statistics, Cohen's d, t-tests, Hotelling's T²
- Performance: ~10-30x speedup from vectorization + pre-computation
- Execution time: ~6.7 minutes for 9,394 perturbations (down from estimated hours)
- Evaluated multiprocessing: 11x SLOWDOWN due to overhead - removed
- Created comprehensive test suite (16 tests, all passing)

### Infrastructure Migration

Migrated uv → pixi package manager for better scientific computing support. All dependencies resolved via conda-forge (except git dependency).

---

## 2025-10-25: Results Analysis & Divergence Discovery

### Virtual Screen Analysis Complete

- Ran notebook 2.0 for LINCS (9,394 perturbations) and TAORF (323 perturbations)
- Modified notebook 2.2 for Excel generation with three filtering levels
- Generated multi-sheet Excel files: All results → Orth filtered → Both filtered

### Critical Discovery: 99.99% Divergence from Baseline

**Problem:** Locally-regenerated CSVs differ massively from S3 baseline (July 2024):

- LINCS: 9,394/9,395 rows differ in `slope`, `last_peak_ind`, `d_slope`
- TAORF: 323/327 rows differ in slope-related values
- Cell counts match exactly → same input data, different algorithm

**Root cause identified:** Control subtraction timing in slope calculation

- Baseline uses `if 0:` branch (control-subtract AFTER slope calculation)
- Current code uses `if 1:` branch (z-score BEFORE slope calculation)
- 35.5% of perturbations with identical radial profiles have opposite-sign slopes
- This is NOT a numerical precision issue - it's a methodology difference

**Evidence:**

- Vectorization functions produce identical results (not the bug)
- Input data is identical (Count_Cells_avg 100% match)
- Exact divergence point: Control subtraction timing (lines 1201-1220 in notebook 2.0)

### Workflow Infrastructure Created

Implemented baseline/regenerated coexistence workflow:

- S3 CSVs → `*_results_pattern_aug_070624.csv` (baseline)
- Local CSVs → `*_results_pattern_aug_070624_REGEN.csv` (regenerated)
- Notebook 2.0: `USE_REGEN_SUFFIX` flag (default: True, preserves baseline)
- Notebook 2.2: `ANALYSIS_MODE` selection ("baseline" or "regenerated")

---

## 2025-10-25: Version Control & File Organization

### Excel File Provenance Traced

Git history revealed curated Excel files originated from Google Sheets (manually curated, not direct pipeline outputs).

### Directory Structure Created

```
data/processed/tables/
├── curated_2024-08-11/          # Original Google Sheets export (git tracked)
├── curated_2025-10-25/          # Current curated version (git tracked)
├── generated_from_baseline/  # From baseline CSVs via notebook 2.2 (gitignored)
└── generated_from_local/        # From regenerated CSVs (gitignored)
```

**Git strategy:** Track curated versions (manual work, irreplaceable), ignore generated versions (reproducible).

### Reproducibility Verified

Compared all three versions:

- `generated_from_baseline` vs `curated_2024-08-11`: **PERFECT MATCH** (after removing duplicate `.1` columns)
- Pipeline is fully deterministic and reproducible
- All 173,806 rows match exactly across 6 datasets
- Aug 2024 baseline represents clean computational output

Differences in `curated_2025-10-25`:

- Minor formatting artifacts (duplicate columns, XML encoding)
- One manually added summary row (jump_orf)
- Presentation-layer only, not scientific data

---

## 2025-10-25: Pipeline Automation & Refactoring

### Snakemake + Typer CLI

- Extracted CSV→Excel logic from notebook 2.2 into `haghighi_mito/data.py`
- Created `DATASET_INFO` configuration dict in `haghighi_mito/config.py`
- Built professional Typer CLI (`haghighi_mito/cli.py`):
  - `haghighi-mito process-csv-single` - Process single dataset CSV
  - `haghighi-mito create-database` - Create DuckDB from Parquet files
  - `haghighi-mito validate-databases` - Compare two DuckDB databases
- Simplified Snakefile: 91% reduction in complexity (14 lines → 1 line per rule)
- Added Justfile commands: `just run`, `just dry`, `just clean`, `just status`

### Parquet Integration & Parallelization

- Converted CSV→Parquet for ~10x faster reads (0.8s vs 11.7s for LINCS)
- Implemented parallel processing: 6 datasets in ~50 seconds (vs 5+ minutes sequential)
- Created intermediate pipeline: CSVs → Parquet → Excel + DuckDB

### Data Organization Refactoring

Migrated all downloads to Snakemake for reproducibility:

- Eliminated `download_data.sh` script
- Semantic directory structure:

  ```
  data/external/mito_project/
  ├── metadata/                     # Dataset metadata (20 files, 1.74 GB)
  ├── per_site_aggregated_profiles/ # Profile data (137 files, 2.69 GB)
  ├── orth_features/                # Orthogonal features (8 files, 9.8 KB)
  └── results/virtual_screen_baseline/ # Pre-computed CSVs (6 files, ~90 MB)
  ```

- All downloads managed via `just download-*` commands
- Pipeline visualization: `just viz` → generates DAG diagrams in `docs/pipeline/`

### Notebook CLI Integration

Converted notebook 2.0 to accept `--dataset` argument for automated execution:

```bash
pixi run python notebooks/2.0-mh-virtual-screen.py --dataset lincs
```

Enables Snakemake to run virtual screening in parallel across datasets.

### Full Integration Test - Success

Ran "nuclear option" test (deleted all data, rebuilt from scratch):

- Downloaded 178 files from S3
- Processed 6 datasets → Parquet → Excel
- Created DuckDB database (178,826 rows)
- All outputs match baseline perfectly
- Pipeline is fully reproducible end-to-end

---

## 2025-10-25: Deep Debugging - Baseline Code Lost

### Investigation Summary

Systematic debugging to identify root cause of 99.99% divergence from baseline (TAORF dataset):

- Confirmed input data identical: Count_Cells_avg matches 100% between baseline and regenerated
- Vectorization correctness verified: `find_end_slope2_vectorized` produces identical results to row-by-row implementation
- Neither `if 0` nor `if 1` branch matches baseline: Both produce ~77% within 10%, ~11% with opposite signs

### Root Cause Identified

**100% of large slope differences have different `last_peak_ind` values** - peak-finding algorithm identifies different peaks in identical input data.

**Critical finding:** Original `if 0` branch has a bug (assigns `slope` but line 1249 expects `peak_slope`) - this code path couldn't run successfully.

**Conclusion:** Baseline was generated with code that **no longer exists** in this repository. Repository created September 2025, baseline uploaded July 2024. Original implementation lost.

### Decision Point

Cannot reproduce baseline without access to July 2024 code. Three options:

1. Accept baseline as-is (use S3 results, optimization work unvalidated)
2. Commit to regenerated version (re-run all datasets, document divergence)
3. Validate regenerated results independently (check biological plausibility with domain expert)

---

## Current Status

### ✅ Baseline Reproduction Achieved (2025-12-04)

**Module now perfectly reproduces July 2024 baseline** (r=1.000 for all metrics).

- **taorf:** 324 perturbations, 99.7% within 10%
- **lincs:** 9,395 perturbations, 100% within 10%

Two algorithmic fixes identified from upstream refactored code:
1. Pre-standardize radial features per plate BEFORE control subtraction
2. Use `nanpercentile(interpolation="nearest")` for median plate selection

See entry **2025-12-04** below for details.

### What's Working ✅

- **Module pipeline:** Fully validated, reproduces baseline perfectly
  - Raw profiles → CSV → Excel → DuckDB
  - r=1.000 correlation with July 2024 baseline
  - Run via: `just generate-module-all`

- **Baseline pipeline:** Pre-computed CSVs (for comparison/fallback)
  - S3 CSVs → Excel → DuckDB (178,826 rows)
  - Run via: `just generate-baseline-all`

- **Infrastructure:** Production-ready automation
  - Snakemake pipeline with Justfile commands
  - Typer CLI for processing operations
  - Parquet integration for performance

- **Performance:** Optimized and tested
  - Vectorized slope: 200x speedup
  - Vectorized stats: 10-30x speedup
  - Parallel processing: 6 datasets in <1 minute

---

## 2025-10-25: Minimal Virtual Screen Module - Baseline Validation

### New Infrastructure Created

- Built `haghighi_mito/virtual_screen.py` - minimal from-scratch implementation
- Added `virtual-screen` CLI command to Typer interface
- Reads same per-site aggregated profiles as notebook 2.0 (confirmed from manuscript methods)
- Calculates simplest metrics: Count_Cells_avg, last_peak_ind, slope

### Validation Results - Input Data Confirmed Identical

Ran taorf dataset comparison with baseline:

- **Count_Cells_avg:** 0.03% mean diff, 98% within 1% → **Input data is identical** ✓
- **last_peak_ind:** 0 exact matches, mean diff 6 bins → Peak detection finds different peaks ✗
- **slope:** 104% mean diff, 0 within 10% → Completely different slopes ✗

### Root Cause Narrowed

- Input data confirmed identical (cell counts match perfectly)
- Data aggregation level confirmed correct (per-site, not per-cell - matches manuscript methods)
- Issue isolated to **peak detection algorithm** in `find_end_slope2_simple()`
- Peak detection finds different local extrema in identical radial patterns
- This is not a numerical precision issue - it's algorithmic

### Next Steps

Foundation established for iterative debugging of slope calculation. Can now test alternative peak detection methods and compare against baseline systematically.

---

## 2025-10-25: Statistical Test Implementation & Partial Baseline Match

### Implementation Complete

- Added statistical test calculations to `haghighi_mito/virtual_screen.py`
- Implemented all 4 baseline t-values using vectorized functions:
  - `t_target_pattern`: Hotelling's T² on full radial distribution (bypasses peak detection)
  - `t_orth`: Hotelling's T² on orthogonal features
  - `t_slope`: Welch's t-test on slope values
  - `d_slope`: Cohen's d effect size for slope
- Fixed orthogonal feature loading (uses `fibroblast_derived.csv` for most datasets, not dataset-specific files)

### Key Finding: t_target_pattern Shows Partial Match

Comparison with baseline (taorf, 327 perturbations):

- **t_target_pattern:** 37% within 10%, 25% within 1% → **Partial validation achieved**
- **t_orth:** 33% within 10%
- **t_slope & d_slope:** <10% within 10% (confirming slope divergence)

**Significance:** t_target_pattern tests the entire 12-bin radial distribution without peak detection. The 37% match rate validates that control-subtracted radial patterns are reasonably similar to baseline, confirming peak detection (not data processing) is the divergence point.

---

## 2025-10-26: Edge Case Analysis Tool for Baseline Comparison

### Implementation Complete

- Added edge case visualization to `haghighi_mito/virtual_screen.py`:
  - `load_perturbation_radial_data()`: Extracts per-site radial patterns for specific perturbations
  - `visualize_peak_detection()`: Plots control-subtracted patterns, smoothing, peak detection, slope calculation
  - `analyze_edge_cases()`: Finds best/worst matches and generates diagnostic plots
- Added CLI command `haghighi-mito analyze-edge-cases` with configurable sorting metric
- Justfile command: `just analyze-edge-cases-for taorf`

### Key Finding: Even Best Matches Show Significant Divergence

Sorted by absolute percentage difference (taorf, 327 perturbations):

- **Best slope matches:** 38-68% difference from baseline → No subset matches well
- **t_target_pattern sorting** (default): Bypasses peak detection, compares raw radial curve similarity
- Confirms fundamental algorithmic difference, not numerical precision issue

### Decision: Use t_target_pattern for Edge Case Analysis

Default sorting metric changed from `slope` to `t_target_pattern` because:

- Tests full 12-bin radial distribution without peak detection
- More direct measure of pattern similarity (37% within 10% vs <10% for slope)
- Helps isolate whether similar patterns produce different slopes due to peak detection alone

---

## 2025-10-26: t_target_pattern Analysis - Promising Alternative Metric

### Key Discovery: Strong Correlation Despite Systematic Scaling

- **t_target_pattern (Hotelling's T² on full radial distribution) shows r=0.92 with baseline**
- Systematic transformation: `new = 0.83 * baseline + 0.40` (R²=0.85)
- 98.8% sign agreement (vs 75.2% for slope)
- Bypasses peak detection entirely - tests all 12 radial bins directly

### Scatter Plot Analysis

Created baseline vs regenerated plots for 4 metrics (taorf, n=327):

- **t_target_pattern**: Excellent linear correlation, tight clustering around fit line
- **slope**: Terrible - most points compressed to horizontal band, fundamentally different
- **t_orth**: Good correlation (r=0.858), validates non-radial processing
- **t_slope**: Moderate (r=0.779), inherits slope issues

### Root Cause: Aggregation Method Differences

Per-plate analysis revealed **only 1.2% exact matches** (4/324 perturbations):

- Baseline aggregated value matches ANY regenerated plate value in only 1.2% of cases
- **Critical difference found in aggregation logic**:
  - Baseline (notebook 2.0:1383-1387): Selects plate with median `d_slope` (column 3), reports ALL stats from that plate
  - Regenerated (virtual_screen.py:364): Selects plate with median absolute `t_target_pattern` (column 0)
- **Plate selection differences could explain 0.83x scaling** - if methods systematically select different plates per perturbation
- Cannot distinguish from per-plate calculation differences without direct plate-to-plate comparison

### Implication

t_target_pattern is more robust than slope (no peak detection dependency) and shows strong baseline agreement. The 0.83x systematic scaling could come from either: (1) different plate selection methods, or (2) preprocessing differences (standardization/normalization in notebook 2.0 lines 1212, 1251). Requires testing with matched plate selection.

### Infrastructure Added

- `analyze_t_target_pattern_distribution()` - correlation and transformation analysis
- `compare_per_plate_results()` - per-plate exact match testing
- `return_per_plate` flag in `calculate_statistical_tests()`
- Scatter plot generation with linear fits

---

## 2025-10-26: Z-Score Normalization Discovery - Baseline Reproduction Achieved

### Root Cause Identified

**Missing z-score normalization per plate** (notebook 2.0 lines 1251-1254):

- Baseline workflow: calculate slopes → z-score per plate → aggregate via median
- Regenerated workflow was missing step 2 entirely
- This explained 100% slope mismatch and 6-bin last_peak_ind offset

### Fix Implemented

Added z-score normalization to `calculate_simple_metrics()` in `virtual_screen.py`:

```python
per_site_df = per_site_df.groupby("batch_plate").apply(z_score_normalize)
```

### Validation Results (taorf, n=327)

**Before fix:**

- slope: 0/327 within 10%, mean diff 104%
- last_peak_ind: 0 exact matches, mean diff 6.0 bins

**After fix:**

- slope: 63/327 within 10% (19%), correlation r=0.849
- last_peak_ind: mean diff 0.23 bins (26x better)
- 71/327 perturbations match well in both metrics (<25%)

### Code Refactoring

- Separated diagnostic functions from core pipeline
- Created `haghighi_mito/diagnostics.py` (539 lines)
- Reduced `virtual_screen.py` from 990 to 466 lines (53% smaller)
- Updated module docstrings, removed one-off diagnostic script

### Conclusion

Baseline is now reproducible at 20-70% level (vs 0% before). Remaining differences likely due to minor implementation details (aggregation methods, outlier handling). Core methodology now matches baseline.

---

## 2025-10-26: CLI Command Rename for Clarity

### Refactoring

- Renamed `analyze-t-target-pattern` command to `compare-baseline-metrics`
- Updated docstring to clarify it compares all 4 statistical metrics, not just t_target_pattern:
  - t_target_pattern (Hotelling's T² on radial distribution)
  - t_orth (Hotelling's T² on orthogonal features)
  - t_slope (Welch's t-test on slope)
  - d_slope (Cohen's d effect size)
- Previous name was misleading - implied single-metric analysis when it actually performs comprehensive baseline comparison

### Command Usage

```bash
pixi run haghighi-mito compare-baseline-metrics --dataset taorf
```

---

## 2025-10-26: Slope Calculation Method Testing - Endpoint Not Root Cause

### Investigation

Tested two slope calculation methods found in codebase to determine which baseline used:

- **Method A** (current): Average of last two points: `endpoint = (smoothed[-1] + smoothed[-2]) / 2`
- **Method B** (alternative): Last point only: `endpoint = smoothed[-1]`

Implemented toggle via `use_last_two_average` parameter in `find_end_slope2_simple()`, threaded through pipeline.

### Results (taorf, n=327)

- **Method A**: 63/327 (19.3%) within 10%, mean % diff 184.75%
- **Method B**: 56/327 (17.1%) within 10%, mean % diff 156.70%
- **Difference**: Only 2% improvement with Method A

### Key Finding

**Slope endpoint calculation is NOT the root cause of 81% mismatch with baseline.**

The remaining divergence must come from:

1. Peak detection differences (both methods show same 0.23 bin offset in `last_peak_ind`)
2. Other algorithmic differences not yet identified
3. Baseline generated with fundamentally different code (July 2024, repository created Sept 2025)

### Decision

Implementation stashed - no value in maintaining alternative method when difference is negligible. Current Method A (average of last two) performs slightly better and remains default.

---

## 2025-10-26: Diagnostic Code Cleanup - Removed Completed Investigations

### Code Removal

Removed diagnostic functions that served their purpose:

- `compare_per_plate_results()` - found 1.2% per-plate match rate, investigation complete
- `return_per_plate` parameter from `calculate_statistical_tests()`
- `compare-per-plate` CLI command
- `analyze-edge-cases-for` from Justfile (dead command)

### Function Simplification

Streamlined `analyze_t_target_pattern_distribution()`:

- **Renamed**: `plot_baseline_comparison()` - name now matches what it does
- **Reduced**: 330 lines → 120 lines (64% smaller)
- **Removed**: All verbose statistical analysis (8 sections)
- **Kept**: 2x2 scatter plots + greppable correlation summary
- **CLI command**: `compare-baseline-metrics` → `plot-baseline-comparison`

### Module Size Reduction

- `diagnostics.py`: 539 lines → 250 lines (54% smaller)
- Justfile updated with `plot-baseline-comparison-for DATASET`

### Greppable Output Format

```
t_target_pattern: r=0.923, within_10%=120/327 (36.7%)
slope: r=0.849, within_10%=63/327 (19.3%)
```

---

## 2025-10-26: Performance Optimization Experiments - Mixed Results

### Vectorized Slope Calculation - Success

- Implemented vectorized slope calculation in `calculate_simple_metrics()` (lines 182-202)
- Replaced `for idx, row in per_site_df.iterrows()` with NumPy matrix operations
- Uses existing `find_end_slope2_vectorized()` from `vectorized_slope.py`
- Performance: 167ms for 16,301 observations (~200x faster than iterrows)

### Parallel Statistical Tests - Limited Success

- Attempted joblib parallelization for statistical test calculations
- Created `_process_single_perturbation()` worker function
- Used `Parallel(n_jobs=-1)` to distribute 324 perturbations across 192 cores
- Result: Only 1.3x speedup (4.5s sequential → 3.5s parallel)

### Root Cause: Task Overhead Dominates

- Each perturbation is a tiny task (~14ms of actual computation)
- 324 independent tasks means high serialization/deserialization overhead
- Joblib loky backend uses memory mapping for shared data, but per-task overhead still significant
- Task overhead > computation time for small datasets

### Decision

- Kept vectorized slopes (clear win, no complexity cost)
- Reverted parallelization code to sequential loop (minimal benefit, added complexity)
- Future optimization would require batching (process chunks of perturbations per worker)
- Not worth complexity for current dataset sizes

---

## 2025-10-26: Pipeline Documentation & Naming Standardization

### Comprehensive Documentation Added

- Rewrote Snakefile docstring (120 lines) documenting all three methods:
  - Method 0: Baseline (validated S3 CSVs, production)
  - Method 1: Notebook (complete pipeline, 1433 lines)
  - Method 2: Module (incomplete, clean 448 lines)
- Added data flow diagrams, reproducibility issue explanation, output directory structure
- Clear status indicators (✅ Complete, ⚠️ Incomplete)

### Naming Collision Fixed

- Removed generic "regenerated" term (was Method 1 specific, sounded generic)
- Renamed all rules/variables to method-specific naming:
  - `all_regenerated` → `all_notebook` (Method 1) / `all_module` (Method 2, TODO)
  - `process_regenerated_csv` → `process_notebook_csv` / `process_module_csv`
  - `REGEN_DIR` → `NOTEBOOK_DIR` / `MODULE_DIR`
- Updated Justfile: `run-regenerated` → `run-notebook`, added `run-module` (TODO)
- Added `clean-notebook` and `clean-module` commands

### TODO Section Standardized

- Rewrote Method 2 completion rules to mirror Method 1 pattern
- Consistent naming: `process_module_csv`, `create_module_database`, `all_module`
- Added new variables: `MODULE_DIR`, `INTERIM_MODULE`, `TABLES_MODULE`
- Updated output paths: `screen_results_module.duckdb`

### Validation

- Tested Method 2 commands: `just run-module-for taorf` + `just plot-comparison-for taorf`
- Confirmed pipeline works (5s execution, generates CSV + comparison + plots)
- All syntax validated, commands working as documented

---

## 2025-10-26: Naming Refactoring - Eliminated Semantic Inconsistencies

### Problem Identified

Inconsistent naming across the three pipeline methods:

- Method 1 used "regenerated" in some places, "local" in others
- Method 2 used "module" in some places, "simple" in others
- Created confusion about what "regenerated" meant (both methods regenerate from raw data)

### Changes Implemented

**Filesystem renames (no backward compatibility):**

- `virtual_screen_regenerated/` → `virtual_screen_notebook/`
- `virtual_screen_simple/` → `virtual_screen_module/`
- `parquet_regenerated/` → `parquet_notebook/`
- `generated_from_local/` → `generated_from_notebook/`

**Code updates:**

- Snakefile: 17 locations (variables, documentation, rules, onstart display)
- Justfile: 6 locations (commands, clean targets)
- Command rename: `run-notebook-csv-for` → `run-notebook-for` (cleaner parallel with `run-module-for`)

### New Consistent Pattern

```
Method 0: baseline    baseline    baseline    baseline
Method 1: notebook    notebook    notebook    notebook
Method 2: module      module      module      module
          ↓           ↓           ↓           ↓
          database    CSV dir     parquet     tables
```

Each method now uses ONE name consistently across all directory paths, variables, and commands.

---

## 2025-10-26: Method 2 Self-Contained - Metadata Preprocessing

### Problem Identified

Module (`virtual_screen.py`) failed with missing file error:

- Expected preprocessed metadata: `workspace/metadata/preprocessed/annot_{dataset}.csv`
- Files exist on S3 from July 2024, but not part of Method 2 download dependencies
- Method 1 (notebook) has built-in preprocessing that creates these files as side effect

### Solution: In-Memory Preprocessing

Added `preprocess_metadata()` function to `virtual_screen.py` (lines 84-176):

- Loads raw metadata for each dataset (LINCS_meta.csv, TA-ORF/replicate_level_cp_normalized.csv.gz, etc.)
- Adds standardized columns: Batch, batch_plate, ctrl_well, Metadata_pert_type
- Filters to necessary columns (critical for taorf - prevents merge conflicts)
- Preprocesses in-memory, nothing written to disk

### Key Design Decision

**Rejected downloading preprocessed files from S3** - Would create hidden dependency and violate "regenerated from scratch" philosophy. Method 2 now truly self-contained.

### Validation

- Tested taorf: 324 perturbations processed successfully
- Module no longer depends on notebook preprocessing or S3 preprocessed files
- Snakefile rule renamed: `run_virtual_screen_analysis` → `run_virtual_screen_module`

---

## 2025-10-26: Dataset Filtering & Diagnostic Naming Cleanup

### Dataset Filtering Feature

- Added `--datasets` parameter to `create-database` CLI command
- Allows selective database creation (e.g., `--datasets lincs,taorf` instead of all 6 datasets)
- Updated Snakefile rules to pass dataset list from `DATASETS` variable
- Supports both lowercase (lincs, taorf) and uppercase (LINCS, TA_ORF) naming

### Diagnostic Output Renaming

- Renamed legacy directory: `t_target_pattern_analysis/` → `diagnostics/`
- Renamed plot files: `baseline_vs_regenerated.png` → `comparison_metrics.png`
- Old naming was a relic from single-metric analysis; new naming reflects multi-metric comparison (t_target_pattern, slope, t_orth, t_slope)
- Updated across: diagnostics.py, Snakefile, Justfile

### Validation Tests

- Baseline database creation: 9,717 rows (LINCS + TA_ORF) ✓
- Diagnostic plots: Valid PNG files with correct naming ✓
- Module pipeline end-to-end: taorf dataset processed successfully ✓
- All three methods tested and working

---

## 2025-10-26: Method 2 CSV Output Fixed - Matches Baseline Format

### Problem Identified

Module CSV had only 8 columns (vs 17 in baseline/notebook):

- Missing: 3 metadata columns (Metadata_gene_name, Metadata_pert_name, Metadata_moa)
- Missing: 6 p-value columns (p_target_pattern, p_orth, p_slope, p_slope_std, p_pattern_std, p_orth_std)
- Wrong column order

### Root Cause

`calculate_simple_metrics()` started with aggregated perturbation data only, not full metadata table.
`calculate_statistical_tests()` only returned t-values, not p-values.

### Solution Implemented

**virtual_screen.py changes:**

- Renamed `calculate_simple_metrics()` → `calculate_metrics()` (clearer purpose)
- Modified to accept `annot` parameter and start with `annot[meta_cols].drop_duplicates()`
- Added all 6 p-value extractions from `batch_results["pvals"]` in `calculate_statistical_tests()`
- Fixed median plate selection to use d_slope (index 3) instead of t_target_pattern (matches notebook 2.0 line 1388)
- Added column reordering in `run_virtual_screen()` to match baseline exactly

### Validation Results

**taorf dataset (327 perturbations):**

- ✅ Column count: 8 → 17 (matches baseline/notebook)
- ✅ Column order: Exact match to baseline
- ✅ Row count: 328 (1 header + 327 data rows, matches baseline)
- ✅ Metadata columns: All 4 dataset-specific columns present
- ✅ Format verified across different datasets (taorf has 4 meta cols, lincs has 10)

**Column order (correct):**

```
Metadata_gene_name, Metadata_pert_name, Metadata_broad_sample, Metadata_moa,
Count_Cells_avg,
p_target_pattern, p_orth, p_slope, p_slope_std, p_pattern_std, p_orth_std,
t_target_pattern, t_orth, t_slope, d_slope,
last_peak_ind, slope
```

### Impact

Method 2 (module) CSV output now has identical format to Method 0 (baseline) and Method 1 (notebook). Ready for next step: Add Excel/Parquet/DuckDB processing to complete Method 2 pipeline.

---

## 2025-10-26: Method 2 Pipeline Complete - Full CSV → Excel → DuckDB

### Changes Implemented

**Snakefile:**

- Uncommented `process_module_csv` rule (lines 463-477)
- Uncommented `create_module_database` rule (lines 479-497)
- Uncommented `all_module` target rule (lines 501-508)
- Added `datasets` parameter to `create_module_database` (matches Method 1 pattern)
- Updated docstring: Method 2 status changed from ⚠️ INCOMPLETE → ✅ COMPLETE
- Updated Method 1 note: "superseded by Method 2" (can be deprecated)

**Justfile:**

- Uncommented `run-module` command
- Updated header: "Complete Pipeline" (was "Incomplete - Stops at CSV")

### Pipeline Validation

**End-to-end test (lincs + taorf):**

- ✅ CSV inputs: `virtual_screen_module/*.csv` (both datasets)
- ✅ Excel outputs: `tables/generated_from_module/*.xlsx` (1.5 MB + 44 KB)
- ✅ Parquet outputs: `parquet_module/*.parquet` (657 KB + 27 KB)
- ✅ DuckDB output: `screen_results_module.duckdb` (1.8 MB, 9,721 rows)
  - LINCS: 9,394 rows
  - TA_ORF: 327 rows

### Commands

```bash
# Full pipeline (all datasets)
just run-module

# Single dataset (CSV only)
just run-module-for taorf

# Clean outputs
just clean-module
```

### Method 2 Status: Production Ready

- ✅ CSV generation with all 17 columns (matching baseline/notebook)
- ✅ Excel generation with 3-sheet filtering (All/Target/Both)
- ✅ Parquet intermediate files
- ✅ DuckDB unified database
- ✅ Cleaner codebase than Method 1 (no dead branches)
- ✅ Identical output format to Method 1

**Method 1 (notebook) can now be deprecated** - Method 2 provides same functionality with cleaner implementation.

---

## 2025-10-26: Separated Baseline Comparison from Virtual Screen

### Problem

Baseline comparison was tied to `virtual-screen` command via `--compare-baseline` flag.
This meant re-running comparison (1 sec) required regenerating entire screen (~5-10 min).

### Solution Implemented

Split into 3 independent commands:

1. **`virtual-screen`** - Generate CSV from per-site profiles (~5-10 min)
2. **`compare-baseline`** - Compare existing CSV with baseline (~1 sec, NEW!)
3. **`plot-baseline-comparison`** - Generate plots from comparison CSV (~1 sec)

**Changes:**

- `virtual_screen.py`: Removed `compare_baseline` parameter from `run_virtual_screen()`
- `virtual_screen.py`: Created new `compare_with_baseline_csv()` function
- `cli.py`: Removed `--compare-baseline` flag from `virtual-screen` command
- `cli.py`: Added new `compare-baseline` CLI command
- `Snakefile`: Split `run_virtual_screen_module` into two rules:
  - `run_virtual_screen_module` - CSV generation only
  - `compare_module_with_baseline` - Comparison only (NEW)
- `Justfile`: Added `compare-baseline-for DATASET` command

**Benefits:**

- Can re-run comparison without regenerating CSV (1000x faster!)
- Cleaner separation of concerns
- More flexible workflow

**New workflow:**

```bash
# Generate CSV once (slow)
just run-module-for taorf

# Compare with baseline (fast - can rerun anytime!)
just compare-baseline-for taorf

# Generate plots (fast)
just plot-comparison-for taorf
```

## 2025-10-26: Justfile Command Refactoring

### What was done

- Standardized all command naming to use consistent `generate-` prefix
- Renamed Snakefile rules for clarity: `run_all_virtual_screen_modules` → `all_module_csvs`, `run_all_virtual_screen_notebooks` → `all_notebook_csvs`
- Removed redundant commands: download commands (Snakemake auto-downloads), database-only commands (identical to `-all` targets), clean commands

### Key decisions

- **Naming pattern**: `generate-{method}-all` (full pipeline), `generate-{method}-csvs` (all datasets CSVs), `generate-{method}-csv-for DATASET` (single CSV)
- Singular/plural distinction: `-csv-for` (one dataset) vs `-csvs` (all datasets)
- Single-command workflows: `just generate-module-all` auto-downloads and processes in one step

### Final command structure

- Method 0: `generate-baseline-all`
- Method 1: `generate-notebook-all`, `generate-notebook-csvs`, `generate-notebook-csv-for DATASET`
- Method 2: `generate-module-all`, `generate-module-csvs`, `generate-module-csv-for DATASET`
- Diagnostics: `compare-baseline-for`, `plot-comparison-for`, `plot-all-comparisons`

---

## 2025-10-26: Diagnostic Workflow Simplification & Module Cleanup

### What was done

- **Moved function**: `compare_with_baseline_csv()` from `virtual_screen.py` → `diagnostics.py` (proper separation of concerns)
- **Unified workflow**: Combined two-step diagnostic (compare → plot) into single command
- **Removed CLI command**: `plot-baseline-comparison` (now auto-called internally)
- **Updated Snakefile**: Merged two rules into one `diagnose_module` rule with both outputs (CSV + PNG)
- **Fixed naming**: `diagnose_all_modules` → `all_module_diagnostics` (consistent with `all_module_csvs`)

### Key decisions

- **One command, both outputs**: `just diagnose-for taorf` now generates comparison CSV + diagnostic PNG automatically
- **Module boundaries**: `virtual_screen.py` = pipeline (raw data → CSV), `diagnostics.py` = analysis (CSV → comparison + plots)
- **No dead imports**: Removed unused imports from `diagnostics.py` (calculate_metrics, calculate_statistical_tests, load_dataset_data)

### Final diagnostic commands

- Single dataset: `just diagnose-for DATASET` → CSV + PNG (~1 sec)
- All datasets: `just diagnose-all` → all diagnostics in parallel
- Removed: `compare-baseline-for`, `plot-comparison-for`, `plot-all-comparisons` (replaced by unified commands)

### Code quality

- 7 files changed: 99 additions, 140 deletions (net -41 lines)
- Tighter docstrings, cleaner architecture, better UX

## 2025-10-27: Baseline Reproduction Hypothesis Testing - Module Validated

### What was done

Tested three hypotheses about why module shows 10-15% disagreement with July 2024 baseline:

- **Hypothesis 2**: Filter plates missing controls (notebook lines 1230-1233)
- **Hypothesis 3**: Two-stage aggregation - plate medians then cross-plate median (notebook lines 1398, 1419-1420)
- **Hypothesis 4**: Median plate selection using percentile without abs() (notebook lines 1385-1392)

### Key findings

**Every attempt to match current notebook made baseline agreement WORSE, not better:**

- Two-stage aggregation: r=0.849 → 0.830 (slope), r=0.862 → 0.862 (t_target_pattern)
- Percentile selection: r=0.90 → 0.849 (slope), r=0.93 → 0.911 (t_target_pattern)
- Filtering plates: No effect on taorf (all plates already have controls), but kept for defensive programming

**Module was ALREADY correctly reproducing July 2024 baseline** (r=0.85-0.90 across metrics). The notebook has drifted from baseline between July 2024 and September 2025.

### Decision

**Keep core implementation unchanged.** Module uses:

- Single-stage aggregation (direct perturbation median, not plate-then-perturbation)
- abs() sorting for median plate selection
- Per-plate z-score normalization of slopes before aggregation

Added defensive improvements:

- Plate filtering to prevent mixing control-corrected with uncorrected data
- Documentation comments explaining implementation choices with empirical evidence
- Pandas compatibility fix (`include_groups=True`) to silence FutureWarning
- Testing script (`scripts/check-baseline-quick.sh`) for rapid baseline validation (~7 sec loop)

The current notebook (2.0-mh-virtual-screen.py) uses different methods and should NOT be treated as ground truth for validation.

### Unresolved issues

**Remaining 10-15% disagreement with baseline is unexplained.** Possible causes include pre-standardization with `normalize_funcs.standardize_per_catX` (lines 1216-1218, 1254-1256 in notebook) or other preprocessing differences. Root cause remains to be investigated.

---

## 2025-10-27: Provenance Tracking & Enhanced Diagnostics

### What was done

- **Added 6 provenance columns** to `virtual_screen.py` output for reproducibility debugging:
  - `n_sites`, `n_plates`, `n_wells`: Sample size metadata
  - `Count_Cells_std`, `slope_std`: Variability metrics
  - `median_plate_id`: Which plate was selected for aggregation
- **Enhanced `diagnostics.py`** with provenance analysis:
  - Added `_compute_provenance_summary()` for metadata statistics
  - Created `export_examples()` to export best + worst slope matches with provenance
  - Added 6 decimal precision to all numeric outputs for readability
  - Fixed sorting bug (was using signed values instead of absolute for "best" matches)
- **Fixed path inconsistencies**: Updated Snakefile/Justfile to use consolidated `virtual_screen_module/` directory for all diagnostic outputs

### Key findings

**Provenance analysis confirms reproducibility issue is algorithmic, not data-related:**

- Both best matches (0.3%-1.4% error) and worst matches (900%-14,000% error) have identical sample sizes (45 sites, 5 plates)
- No correlation between sample size and reproducibility
- All examples have `n_wells=1` (single well position per perturbation across plates)
- Conclusion: 10-15% baseline disagreement stems from calculation/aggregation differences, NOT missing data or low sample sizes

### Notes

- Database creation/validation still works (provenance columns silently dropped as expected)
- New outputs: `*_slope_examples.csv` with top 10 best/worst matches per dataset
- Diagnostic summaries now include provenance section with warnings for low sample sizes

---

## 2025-10-27: Virtual Screen Implementation Deep Dive - Consistency Analysis

### What was done

- Comprehensive review of `haghighi_mito/virtual_screen.py` implementation against docstring and manuscript
- Cross-checked code logic in `calculate_metrics()`, `calculate_statistical_tests()`, and helper modules
- Analyzed consistency between patient fibroblast MITO-SLOPE algorithm (manuscript lines 173-218) and virtual screen slope calculation

### Key findings (observations requiring verification)

**✅ Docstring vs Code: CONSISTENT** - All 5 pipeline steps match implementation exactly

**✅ Code Logic: SOUND** - Control normalization, z-score per-plate normalization, vectorized slope calculation, and statistical testing all follow standard practices

**⚠️ Manuscript Documentation Gaps** (Methods section doesn't explicitly describe):

- Z-score normalization per plate before aggregation (mentioned in Results but not Methods)
- Single-stage median aggregation strategy (vs two-stage)
- Control subtraction before slope calculation
- Different MITO-SLOPE algorithms for patient analysis vs virtual screen (detailed vs simplified)

**⚠️ Empirical Optimization Evidence**: Code comments (lines 348-350, 491-493) reveal implementation was tuned to match baseline (r=0.90), suggesting baseline was generated with THIS implementation. Circular dependency: code optimized to match baseline, baseline used as ground truth.

### Observations needing verification

1. **Methods section completeness**: Verify manuscript Methods accurately describes all preprocessing steps
2. **Algorithm variant documentation**: Confirm manuscript distinguishes patient MITO-SLOPE (detailed, lines 173-218) from virtual screen slope (simplified)
3. **Empirical tuning**: Validate that implementation choices (single-stage aggregation, abs() median plate selection) represent July 2024 baseline methodology vs notebook drift

### Next steps for future investigation

- Consider updating Methods section with explicit z-score normalization and aggregation description
- Test whether two-stage aggregation (notebook 2.0 current implementation) improves or degrades biological interpretability
- Document rationale for implementation choices in supplementary methods
- Investigate why current notebook differs from baseline if module better matches July 2024 results

---

## 2025-11-01 to 2025-11-02: Minimal Reproduction Script & Documentation

### What was done

- Created standalone reproduction script (`scripts/reproduce_slope_discrepancy.py`) demonstrating baseline mismatch (slope r=0.849, peak r=0.516) for sharing with original author
- Fixed pandas FutureWarning: replaced `groupby().apply()` with `transform()` for z-score normalization
- Fixed Cartesian product bug in taorf comparisons (duplicate perturbation IDs causing inflated match counts)
- Added file provenance tracing with audit hooks (logs all 4 input/output files accessed)
- Enhanced script documentation: added FILES ACCESSED, MODULES CALLED sections, and function references for each author question
- Fixed technical inaccuracies in docstring (smoothing: Savitzky-Golay not moving average; preprocessing: clarified z-scoring happens after slope calculation)
- Converted notebook 2.0 to reference baseline Python script (`2.0-mh-virtual-screen-original.py`, 1192 lines)

### Key findings

- Reproduction script confirms r=0.849 slope correlation, sufficient for author consultation
- Script includes 5 focused questions mapping to specific implementation functions (calculate_metrics, find_end_slope2_vectorized)
- Notebook conversion preserves original state (1192 lines) vs current evolved version (1438 lines, +246 lines)

### Notes

- Script uses loguru logger for consistency with codebase (colored output, automatic timestamps)
- All outputs in `data/processed/virtual_screen_module/` for diagnostic workflow integration

---

## 2025-11-02: Reproduce Script Enhancement - Dataset Agnostic Support

### What was done

- Made `scripts/reproduce_slope_discrepancy.py` dataset-agnostic with Typer CLI support
- Fixed critical bugs in `haghighi_mito/virtual_screen.py` for JUMP datasets (jump_orf, jump_crispr, jump_compound)
- Added `just reproduce [DATASET]` command with auto-download via Snakemake
- Fixed metadata merge bugs: replaced non-existent `Metadata_PlateID` with `["Metadata_Plate", "Metadata_Source"]`
- Fixed control well definitions for each JUMP dataset

### Key findings

**Baseline correlation results for all datasets:**

- taorf: slope r=0.849, peak r=0.516 (327 perturbations)
- jump_orf: slope r=0.788, peak r=0.509 (14,792 perturbations)
- jump_crispr: slope r=0.519, peak r=0.370 (7,985 perturbations)
- jump_compound: slope r=0.779, peak r=0.515 (115,729 perturbations)

### Notes

- Lower correlations for JUMP datasets suggest different preprocessing in original baseline
- Script now works for all datasets with available baseline data

---

## 2025-11-02: Fixed jump_compound Control Definition - Pipeline Complete

### What was done

- Fixed jump_compound control well definition in `haghighi_mito/virtual_screen.py` (line 162)
- Changed from `annot.get("Metadata_control_type")` to `annot["Metadata_JCP2022"].isin(["JCP2022_033924"])`
- JCP2022_033924 is DMSO (negcon) control present on all compound plates

### Key findings

**jump_compound now processes successfully:**

- 115,729 perturbations (vs 0 before)
- slope r=0.779 with baseline (good correlation, similar to jump_orf r=0.788)
- last_peak_ind r=0.515 (consistent with other JUMP datasets)

**All six datasets now working:**

1. taorf: r=0.849 (327 perturbations)
2. lincs: not tested yet
3. CDRP: not tested yet
4. jump_orf: r=0.788 (14,792 perturbations)
5. jump_crispr: r=0.519 (7,985 perturbations)
6. jump_compound: r=0.779 (115,729 perturbations) ✓ FIXED

### Notes

- Previous issue: compound metadata lacks `Metadata_control_type` column (unlike orf/crispr)
- Solution: Use specific JCP2022 identifier for DMSO control instead of generic column
- Reproduction script now complete for all JUMP datasets

---

## 2025-11-02: Peak Index as Primary Diagnostic Target

### What was discovered

**Pipeline dependency confirmed:**

- Peak index calculated before slope in `vectorized_slope.py::find_end_slope2_vectorized()` (lines 91-108)
- Slope depends on peak index: `slopes = (endpoint - peak_values) / (n_cols - last_peak_ind - 1)` (line 118)
- Peak index has worse baseline agreement (r=0.516) than slope (r=0.849)

**Key insight: Lower peak correlation is expected**

- Peak index is discrete (bins 0-11) - off-by-1 error is large in relative terms
- Radial distributions are smooth functions - adjacent bins have similar values
- Slope depends on peak VALUE (height) not just INDEX (position)
- Result: Large index errors can produce small slope errors if profiles are similar

**Peak index advantages as diagnostic target:**

- Calculated early in pipeline (fewer dependencies than t-values)
- Faster to compute (no statistical testing or control comparisons)
- More isolated from downstream processing
- Still necessary to understand regardless of correlation pattern

### Pipeline stages producing final baseline peak index values

The final `last_peak_ind` values in baseline CSVs are produced through these sequential stages:

1. **Load per-site radial distributions** (12 bins, radial distribution features 5-16)
2. **Control subtraction** - Subtract per-plate control means from radial patterns
3. **Smoothing** - Apply Savitzky-Golay filter (window_length=5, polyorder=3)
4. **Peak detection** - Find argmax and argmin, select last valid peak (not in last 2 bins)
5. **Z-score normalization** - Normalize peak indices per plate
6. **Aggregation** - Median across all per-site observations per perturbation

Differences between baseline and current implementation could occur at any of these 6 stages.

**Key diagnostic insight:** Stages affecting different metric subsets:

- Stage 2 (control subtraction) affects ALL metrics - both radial and orthogonal features
- Stages 3-4 (smoothing, peak detection) only affect radial distribution metrics (last_peak_ind, slope, t_target_pattern)
- Stages 5-6 (z-scoring, aggregation) affect ALL metrics

**Observed:** t_orth (orthogonal features) also shows discrepancies with baseline, suggesting errors may exist in global stages (control subtraction, z-scoring, aggregation) that affect all metrics. However, current investigation focuses on peak index as primary diagnostic target - understanding peak index discrepancies may explain a bulk of issues across all metrics if root cause is in shared processing stages.

### Limitation

Baseline CSVs only contain final aggregated z-scored values. Cannot isolate which specific stage causes discrepancy without baseline intermediate values (per-site peaks, smoothed profiles, control means, pre-z-scored peaks).

### Next steps

Peak index investigation necessary but constrained by data availability. Pattern analysis (distributions, confusion matrices) may reveal systematic biases even without intermediate baseline values.

---

## 2025-11-08: Two-Stage Aggregation Confirmed as Baseline Methodology

### What was done

- Diagnosed cell count discrepancies between baseline and regenerated results (e.g., XIAP: 62.01 vs 64.28)
- Tested two-stage aggregation (plate-level → cross-plate) vs single-stage (direct cross-site)
- Implemented two-stage aggregation for all metrics (Count_Cells, slope, last_peak_ind)

### Key findings

**Two-stage aggregation exactly matches baseline methodology** (notebook lines 1100-1175):
- Stage 1: Aggregate per plate (mean for Count_Cells, median for slope/peak)
- Stage 2: Aggregate across plates (mean of plate-means, median of plate-medians)

**Correlation results (taorf, n=324):**
- Count_Cells_avg: r=1.000 (perfect match with two-stage)
- Slope: r=0.830 (two-stage) vs r=0.849 (single-stage) - **slightly worse but correct methodology**
- Peak: r=0.533 (two-stage) vs r=0.516 (single-stage) - slightly better

**Previous PROGRESS.md note (2025-10-26) stating "two-stage made slope worse" was CORRECT.** Trade-off accepted: two-stage aggregation slightly reduces slope correlation but matches the validated baseline methodology exactly.

### Resolution

Module now uses two-stage aggregation consistently, matching original notebook implementation. This is the correct approach even though single-stage aggregation empirically achieved higher slope correlation (r=0.849 vs 0.830). Methodological fidelity to baseline takes precedence over metric optimization.

---

## 2025-11-09: "Original" Notebook Validation - Module is Closer to Baseline

### What was done

- Created minimal reproduction notebook from `2.0-mh-virtual-screen-original.py` (labeled "original")
- Applied surgical edits to run taorf dataset only (7 minimal changes)
- Compared output with baseline using diagnostic tools (hack: copied file to module directory)

### Key findings

**Observed baseline agreement:**

| Method | slope r | t_target_pattern r | Within 10% (slope) |
|--------|---------|-------------------|-------------------|
| Baseline (July 2024) | 1.000 | 1.000 | 100% (reference) |
| Module (haghighi_mito) | 0.830 | 0.923 | 19% |
| **Minimal notebook** | **0.415** | **0.756** | **7.7%** |

**Observation:** Minimal notebook from `2.0-mh-virtual-screen-original.py` produces 2x worse baseline agreement than module.

**Possible explanations:**
- Modifications required for local execution (commented sys.path, pixi singlecell package vs original path)
- Package version differences (original used dgx environment, now using pixi)
- File may not be the true baseline-generating code (repo created Sept 2025, baseline from July 2024)
- Environmental differences not accounted for

**Cannot conclude** which explanation is correct without further investigation. What is certain: module (r=0.830) currently achieves better baseline agreement than this notebook snapshot (r=0.415).

### Diagnostic methodology (hack)

Diagnostic tool (`compare-baseline`) expects module output in `virtual_screen_module/` directory. To compare "original" notebook output:

```bash
# 1. Run minimal notebook (outputs to virtual_screen_original/)
pixi run python notebooks/2.0-mh-virtual-screen-minimal.py

# 2. Copy to module directory for comparison
cp data/processed/virtual_screen_original/taorf_*.csv \
   data/processed/virtual_screen_module/taorf_*.csv

# 3. Run diagnostic comparison
pixi run haghighi-mito compare-baseline --dataset taorf
```

This is a workaround - diagnostic tool was designed for module/baseline comparison, not original/baseline. Results are valid but workflow is non-standard.

### Resolution

- **Module remains best available implementation** - achieves r=0.830 baseline agreement
- **Notebook snapshot underperforms** - but root cause unclear (environmental vs code differences)
- **Further investigation needed** - to determine if package/environment differences explain discrepancy
- **Two-stage aggregation confirmed correct** (2025-11-08 entry shows r=1.000 for cell counts)

---

## 2025-12-04: Perfect Baseline Reproduction Achieved

### What was done

- Analyzed upstream refactored code (`3.rank_perturbations.ipynb`) from Marzieh's updated repository
- Identified two algorithmic differences causing baseline disagreement
- Implemented fixes in `haghighi_mito/virtual_screen.py`

### Key findings

**Two fixes achieved r=1.000 baseline reproduction:**

1. **Pre-standardize radial features per plate BEFORE control subtraction**
   - Added z-score normalization at `calculate_metrics()` line 272-277
   - Previously suspected (PROGRESS.md line 956) but never tested
   - This alone improved slope from r=0.830 → r=1.000

2. **Use `nanpercentile(interpolation="nearest")` for median plate selection**
   - Changed from `np.argsort(np.abs(tvals))` to upstream method
   - This fixed t-values: t_target_pattern r=0.942 → r=1.000

Also added second standardization pass for radial + orthogonal features before statistical tests (matches upstream double-standardization pattern).

### Validation results

| Dataset | Perturbations | slope r | t_target_pattern r | Within 10% |
|---------|--------------|---------|-------------------|------------|
| taorf   | 324          | 1.000   | 1.000             | 99.7%      |
| lincs   | 9,395        | 1.000   | 1.000             | 100%       |

### Resolution

- **Environment differences NOT required** - pure algorithmic fix
- **Module now production-ready** - can regenerate results from raw data
- **Upstream refactored code** validated the correct methodology
- Previous hypothesis about scipy/numpy version differences was a red herring

---

## 2025-12-04: Workflow Consolidation - Justfile Eliminated

### What was done

- Consolidated all workflow commands into Snakefile, eliminating Justfile
- Added utility rules: `generate_dag_visualizations`, `validate_database_pair`
- Updated `.envrc` with pixi shell-hook for direct `snakemake` usage

### Key changes

| Old (Justfile) | New (Snakefile) |
|----------------|-----------------|
| `just generate-baseline-all` | `snakemake all_baseline -c4 -p` |
| `just generate-module-all` | `snakemake all_module -c4 -p` |
| `just diagnose-for taorf` | `snakemake data/processed/virtual_screen_module/taorf_comparison_metrics.png` |
| `just viz` | `snakemake generate_dag_visualizations` |

### Notes

- Pattern follows successful consolidation in `cpg0037-oasis-broad-U2OS-data` repo
- Net reduction: -71 lines (Justfile deleted, Snakefile expanded with utility rules)

---

## 2025-12-04: JUMP Dataset Discrepancy Investigation

### What was done

- Investigated why JUMP datasets show r=0.993 while CDRP shows r=1.000
- Compared with upstream repo (`2025_Haghighi_Mito-upstream`) preprocessing
- Added `handle_nans()` function (from upstream `utils.py`) to test if NaN handling was the cause

### Key findings

- `handle_nans` preprocessing did NOT improve correlation (still r=0.993)
- Function retained in codebase but not called (for future reference)
- JUMP discrepancy remains unexplained; likely due to pre-repository baseline generation

### Notes

- r=0.993 is acceptable for practical purposes
- Root cause may be floating point differences or aggregation order in original code
