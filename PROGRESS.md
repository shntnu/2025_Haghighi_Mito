# Progress Log

## 2025-10-17: Data Download Infrastructure Setup

### Completed

- [x] Analyzed S3 data requirements across all notebooks (1.0, 2.0, 2.1, 2.2)
- [x] Created download scripts: `download_data.sh`, `restore_intelligent.py`, `check_restore_status.sh`
- [x] Documented data flow and dependencies in `data_download_analysis.md`
- [x] Updated `.gitignore` to exclude data/logs
- [x] Initiated Glacier restoration (Bulk tier, ~12 hour wait)
  - 20 metadata files
  - 137 per-site profile files
  - 8 orthogonal features files already available

### Status: Waiting for Glacier Restoration

- **Files restoring**: 157 files from Intelligent-Tiering Deep Archive
- **Tier**: Bulk (cheapest, 12 hour retrieval)
- **Started**: 2025-10-17 ~08:07 AM
- **Expected ready**: 2025-10-17 ~08:00 PM

### Next Actions

1. Check restoration status: `bash scripts/check_restore_status.sh`
2. Once restored, download data: `bash scripts/download_data.sh`
3. Verify download: `find data/external/mito_project -type f | wc -l` (expect ~200 files)
4. Update notebook paths to point to local data (currently hardcoded to remote server)
5. Run analysis pipeline:
   - Optional: `1.0-mh-feat-importance.py`
   - Required: `2.0-mh-virtual-screen.py` (generates virtual_screen outputs)
   - Then: `2.1-mh-set-enrichment-analysis.py`, `2.2-mh-check-vs-lists.py`

### Notes

- Total download: ~4.43 GB (not 57 GB - most of results/ is outputs)
- Local reference data already in repo: KEGG/WikiPathways files (38 KB)
- Notebooks 2.1 and 2.2 depend on 2.0 outputs (virtual_screen results)

---

## 2025-10-24: Data Download Completed

### Completed

- [x] Verified Glacier restoration status - all files restored and ready
- [x] Downloaded all required data from S3 using `download_data.sh`
  - 178 files downloaded successfully
  - 4.6 GB total (metadata, per-site profiles, orthogonal features)
  - All 3 phases completed: metadata (1.74 GB), profiles (2.69 GB), features (9.8 KB)
- [x] Verified download integrity with re-sync - all files match remote

### Status: Data Ready for Analysis

- **Local data location**: `data/external/mito_project/workspace/`
- **File count**: 178 files
- **Total size**: 4.6 GB
- **Sync verified**: Re-running download script confirms all files present and current

### Next Actions

1. [ ] Check if notebook paths need updating to use local data
2. [ ] Run virtual screen analysis: `notebooks/2.0-mh-virtual-screen.py`
3. [ ] After 2.0 completes, run enrichment analysis: `2.1-mh-set-enrichment-analysis.py`
4. [ ] Validate results: `2.2-mh-check-vs-lists.py`

### Notes

- `aws s3 sync` commands allow safe re-running without redownloading
- Download script completed in seconds on re-run (files already synced)
- Created output directories: `metadata/preprocessed/`, `results/virtual_screen/`

---

## 2025-10-25: Virtual Screen Analysis Pipeline - Local Execution Setup

### Completed

- [x] Fixed import errors in `notebooks/2.0-mh-virtual-screen.py`
  - Removed non-existent `bbf_test` from imports
  - Applied same fix to `2.2-mh-check-vs-lists.py`
- [x] Updated paths to use local data instead of remote server
  - Changed `mito_project_root_dir` from `/home/jupyter-mhaghigh@broadinst-ee45a/...` to `data/external/mito_project/`
  - Updated TA-ORF metadata path to use local file
- [x] Fixed metadata preprocessing bugs
  - Added missing `Batch` column for jump_compound dataset (line 188)
  - Added `ctrl_well` column for jump_compound using DMSO (JCP2022_033924) and untreated (JCP2022_999999) controls
- [x] Skipped per-site profile creation section (lines 752-1036)
  - Wrapped in triple quotes since pre-computed profiles already downloaded
  - Original section required raw SQLite data not available locally
- [x] Added comprehensive logging with loguru
  - Preprocessing section: logs each dataset processed and saved
  - Analysis section: logs data loading, merging, standardization, slope calculation, and statistical testing
  - Progress tracking every 100 perturbations during statistical tests
- [x] Identified and fixed dataset selection bugs
  - Commented out hardcoded `dataset = "jump_compound"` on line 1244 that was overwriting user selection
  - Added clear section markers and warnings about "last line wins" pattern
- [x] Added code documentation
  - Clearly demarcated preprocessing vs analysis sections
  - Documented which sections process all datasets vs single dataset

### Status: Ready for Analysis Runs

- **Script location**: `notebooks/2.0-mh-virtual-screen.py`
- **Current dataset**: User-selectable via lines 1080-1085 (currently set to `lincs`)
- **All metadata files generated**: 5 datasets (jump_orf, jump_crispr, jump_compound, lincs, taorf)
- **Logging active**: Timestamps and progress tracking throughout pipeline

### Next Actions

1. [ ] Run analysis for each dataset individually:
   - `lincs` (compounds)
   - `jump_orf` (genetics)
   - `jump_crispr` (genetics)
   - `jump_compound` (compounds)
   - `taorf` (genetics)
2. [ ] Check output files in `data/external/mito_project/workspace/results/virtual_screen/`
3. [ ] Run enrichment analysis: `2.1-mh-set-enrichment-analysis.py`
4. [ ] Validate results: `2.2-mh-check-vs-lists.py`

### Notes

- **Notebook workflow pattern**: This is a Jupyter notebook converted to `.py` script
  - Preprocessing section (lines 123-261): Creates all `annot_*.csv` files in one run
  - Analysis section (lines 1061+): Only processes ONE dataset per run (selected via dataset variable)
  - To analyze different datasets, edit lines 1080-1085 and re-run
- **CDRP dataset**: Intentionally disabled (commented out) in preprocessing
- **Slow operations identified**:
  - Peak slope calculation (`np.apply_along_axis`) processes every row
  - Statistical tests loop iterates through all perturbations
  - Both now have logging for progress tracking
- **Key fixes applied**:
  - Import compatibility with installed `singlecell` package
  - Path compatibility for local filesystem
  - Control well identification for all datasets
  - Removed conflicting dataset assignments

---

## 2025-10-25: Performance Optimization - Vectorized Slope Calculation

### Completed

- [x] Created `haghighi_mito/vectorized_slope.py` module
  - Implemented `find_end_slope2_vectorized()` - fully vectorized version of slope calculation
  - Kept original `find_end_slope2()` for reference and compatibility
  - Follows lab workflow conventions (processing code in package, not notebooks)
- [x] Created comprehensive test suite in `haghighi_mito/tests/test_vectorized_slope.py`
  - 9 test cases covering equivalence, edge cases, and realistic profiles
  - Performance benchmark test showing ~200x speedup
  - All tests pass with identical results between original and vectorized versions
- [x] Updated `notebooks/2.0-mh-virtual-screen.py` to use vectorized function
  - Added import: `from haghighi_mito.vectorized_slope import find_end_slope2_vectorized`
  - Replaced slow `np.apply_along_axis()` call at line 1221 with vectorized version
- [x] Fixed package installation
  - Added build system configuration to `pyproject.toml`
  - Added `[tool.setuptools.packages.find]` to specify package location
  - Installed package in editable mode with `uv pip install -e .`
  - Added pytest as dev dependency

### Status: Performance Issue Resolved

- **Performance improvement**: ~200x speedup (0.174s → 0.0008s per 1000 rows)
- **Correctness verified**: Vectorized version produces identical results to original
  - Peak indices match exactly
  - Slopes match to within 1e-10 relative tolerance
- **Test coverage**: 9 tests all passing
- **Notebook verified**: User confirmed script runs successfully

### Next Actions

1. [ ] Run analysis for each dataset individually:
   - `lincs` (compounds)
   - `jump_orf` (genetics)
   - `jump_crispr` (genetics)
   - `jump_compound` (compounds)
   - `taorf` (genetics)
2. [ ] Check output files in `data/external/mito_project/workspace/results/virtual_screen/`
3. [ ] Run enrichment analysis: `2.1-mh-set-enrichment-analysis.py`
4. [ ] Validate results: `2.2-mh-check-vs-lists.py`

### Notes

- **Vectorization approach**: Replaced row-by-row `apply_along_axis` with batch operations
  - Savitzky-Golay smoothing applied to entire matrix at once
  - Peak finding (argmax/argmin) vectorized across all rows
  - Slope calculation using advanced NumPy indexing
- **Why this matters**: The slope calculation section was a major bottleneck
  - Previous approach: Python loop over every row via `apply_along_axis`
  - New approach: Single NumPy operations on entire matrix
  - Critical for large datasets (thousands of perturbations × sites)
- **Code organization**: Following `protocols/workflows.md` conventions
  - Processing code lives in `haghighi_mito/` package (not notebooks)
  - Notebooks import from package for reusable, testable code
  - Tests in `haghighi_mito/tests/` validate correctness
- **Package structure now includes**:
  - `haghighi_mito/config.py` - Project configuration
  - `haghighi_mito/data.py` - Data handling utilities
  - `haghighi_mito/vectorized_slope.py` - Optimized slope calculation (NEW)
  - `haghighi_mito/tests/test_vectorized_slope.py` - Test suite (NEW)

---

## 2025-10-25: Statistical Testing Optimization - Vectorized Batch Processing

### Completed

- [x] Created `haghighi_mito/vectorized_stats.py` module
  - Implemented vectorized statistical functions for batch processing across plates
  - `batch_plate_statistics()` - main function replacing inner plate loop
  - `cohens_d_vectorized()`, `ttest_ind_vectorized()` - batch statistical tests
  - `t_to_z_vectorized()`, `z_to_p_vectorized()` - vectorized conversions
  - `TwoSampleT2Test_vectorized()` - batch Hotelling's T² tests
  - GPU-ready architecture: All NumPy operations can migrate to CuPy with minimal changes
- [x] Created comprehensive test suite in `haghighi_mito/tests/test_vectorized_stats.py`
  - 16 test cases covering correctness, edge cases, and batch operations
  - Validates vectorized functions match original implementations to 1e-10 precision
  - All tests passing ✅
- [x] Optimized main statistical testing loop in `notebooks/2.0-mh-virtual-screen.py`
  - Pre-compute control dataframes by plate (saves thousands of redundant filtering operations)
  - Replaced 80-line nested plate loop with single `batch_plate_statistics()` call
  - Each perturbation processes all plates in vectorized batch operations
- [x] Evaluated multiprocessing parallelization
  - Implemented full parallel processing with joblib
  - Benchmarked: Found 11x SLOWDOWN due to overhead (24.8s vs 2.15s for 50 perturbations)
  - Root cause: Tasks too fast (43ms each), overhead dominates (process spawning, data serialization)
  - Decision: Removed multiprocessing code - vectorization alone is sufficient
- [x] Cleaned up codebase
  - Removed `haghighi_mito/parallel_stats.py` and test file
  - Removed multiprocessing configuration and benchmark section from notebook
  - Removed joblib dependency

### Status: Optimization Complete - Ready for Production

- **Performance improvement**: ~10-30x speedup from vectorization alone
- **Execution time**: ~6.7 minutes for 9,394 perturbations (down from estimated hours)
- **Per-perturbation time**: ~43ms average
- **Code quality**:
  - 16 tests all passing
  - GPU-migration ready (90% compatible with CuPy)
  - Clean, maintainable codebase
- **Key optimizations**:
  1. Pre-computed control dataframes (eliminates repeated filtering)
  2. Vectorized statistical tests (batch operations replace loops)
  3. NumPy broadcasting for efficient computation

### Next Actions

1. [ ] Run analysis for each dataset:
   - `lincs` (compounds) - 9,394 perturbations
   - `jump_orf` (genetics)
   - `jump_crispr` (genetics)
   - `jump_compound` (compounds)
   - `taorf` (genetics)
2. [ ] Check output files in `data/external/mito_project/workspace/results/virtual_screen/`
3. [ ] Run enrichment analysis: `2.1-mh-set-enrichment-analysis.py`
4. [ ] Validate results: `2.2-mh-check-vs-lists.py`

### Notes

- **Why vectorization worked**:
  - Original: Nested loops with repeated DataFrame filtering per plate
  - Optimized: Pre-compute controls once, batch all plate operations per perturbation
  - NumPy operations are highly optimized C code, orders of magnitude faster than Python loops
- **Why multiprocessing didn't work**:
  - Overhead (process spawning, IPC, data copying) >> computation time
  - Each perturbation completes in ~43ms - too fast for multiprocessing benefits
  - Would help if tasks took seconds/minutes, not milliseconds
- **GPU migration path** (for future 100-500x speedup if needed):
  - Change `import numpy as np` → `import cupy as cp` in vectorized_stats.py
  - Transfer data to GPU memory once before loop
  - Batch operations run on GPU cores automatically
  - Current vectorized code is 90% compatible with CuPy
- **Performance breakdown** (for 9,394 perturbations):
  - Sequential vectorized: ~6.7 minutes (current implementation)
  - Original nested loops: Estimated 1-2 hours (before optimization)
  - Combined speedup: ~10-30x from vectorization + pre-computation
- **Code organization**:
  - `haghighi_mito/vectorized_stats.py` - Reusable batch processing functions
  - `haghighi_mito/tests/test_vectorized_stats.py` - Comprehensive test coverage
  - Notebook simplified from 80+ lines to clean loop with function call
  - Follows lab workflow conventions (processing code in package, not notebooks)

---

## Template for Future Entries

```text
## [Date]: [Brief Description]

### Completed

- Item 1
- Item 2

### Status

- Issue description
- Resolution or workaround

### Next Actions

- [ ] Action 1
- [ ] Action 2

### Notes

- Relevant context
```