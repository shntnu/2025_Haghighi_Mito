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