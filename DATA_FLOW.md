# Data Flow Documentation: Mitochondrial Morphology Virtual Screening Pipeline

## Overview

This document traces the complete data flow through the mitochondrial morphology analysis pipeline, from raw patient fibroblast images and public Cell Painting datasets to final virtual screening results identifying chemical and genetic modulators of mitochondrial localization.

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     PATIENT FIBROBLASTS (1.0)                   │
│  Goal: Identify disease-associated mitochondrial phenotype      │
└────────────┬────────────────────────────────────────────────────┘
             │
             ├─→ Feature Discovery → MITO-SLOPE metric
             ├─→ Machine Learning Classification → ROC-AUC
             └─→ Statistical Validation → Figures 3, 4
                          ↓
┌─────────────────────────────────────────────────────────────────┐
│                  PUBLIC DATASETS SCREENING (2.0)                │
│  Goal: Find perturbations affecting MITO-SLOPE in Cell Painting │
└────────────┬────────────────────────────────────────────────────┘
             │
             ├─→ Step 1: Metadata Preprocessing
             ├─→ Step 2: Identify Orthogonal Features
             ├─→ Step 3: Per-Site Aggregation
             └─→ Step 4: Statistical Testing → Results CSVs
                          ↓
┌─────────────────────────────────────────────────────────────────┐
│                   ENRICHMENT ANALYSIS (2.1)                     │
│  Goal: Identify enriched pathways/gene sets in hits             │
└────────────┬────────────────────────────────────────────────────┘
             │
             └─→ GSEA on GO, KEGG, WikiPathways, MOA, etc.
                          ↓
┌─────────────────────────────────────────────────────────────────┐
│                 FINAL HIT SELECTION (2.2)                       │
│  Goal: Apply statistical filters and generate final tables      │
└────────────┬────────────────────────────────────────────────────┘
             │
             └─→ BH correction → Filtered XLSX → Manuscript Tables
```

---

## Phase 0: Raw Data Sources

### Patient Fibroblast Data (Remote S3)

```
Location: ~/bucket/projects/2016_08_01_RadialMitochondriaDistribution_donna/

workspace/
├── backend/Mito_Morphology_input/
│   └── {subject}/
│       └── {subject}.sqlite           # Single-cell CellProfiler features
│
├── singleCellData/
│   ├── single_cell_with_annot.pkl     # Pre-existing single-cell profiles
│   └── single_cell_with_annot_allFeatures.pkl  # Expanded feature set
│
└── metadata/
    └── patient_labels_updatedSept302025.xlsx   # Disease labels (168 patients)
```

**Patient Categories:**
- BP (Bipolar Disorder): 50 patients
- SZ (Schizophrenia): 39 patients
- SZA (Schizoaffective): 11 patients
- MDD (Major Depression): 21 patients
- Control (Healthy): 47 patients

### Public Cell Painting Datasets (Remote)

```
1. LINCS (cpg0004)
   Location: ~/bucket/projects/2015_10_05_DrugRepurposing.../workspace/backend/
   Content: 1,571 compounds × 6 doses in A549 cells

2. CDRP (cpg0012)
   Location: ~/gallery/cpg0012-wawer-bioactivecompoundprofiling/broad/workspace/backend/
   Content: 30,430 compounds in U2OS cells

3. JUMP-ORF (cpg0016-jump/source_4)
   Location: ~/gallery/cpg0016-jump/source_4/workspace/backend/
   Content: 12,609 overexpressed genes in U2OS cells

4. JUMP-CRISPR (cpg0016-jump/source_13)
   Location: ~/gallery/cpg0016-jump/source_13/workspace/backend/
   Content: 7,975 CRISPR knockouts in U2OS cells

5. JUMP-Compound (cpg0016-jump/source_13)
   Location: ~/gallery/cpg0016-jump/source_13/workspace/backend/
   Content: 116,000+ compounds in U2OS cells

6. TA-ORF (cpg0017)
   Location: ~/gallery/cpg0017-rohban-pathways/broad/workspace/backend/
   Content: 323 overexpressed genes in U2OS cells
```

**Data Format:** SQLite databases, CSV files, or Parquet files at path:
```
{backend}/{batch}/{plate}/{plate}.{sqlite|csv|parquet}
```

---

## Notebook 1.0: Patient Fibroblast Feature Discovery

**File:** `notebooks/1.0-mh-feat-importance.py`

### Input Data

```python
INPUT FILES:
├── single_cell_with_annot.pkl                    # Pre-existing profiles
├── single_cell_with_annot_allFeatures.pkl        # Expanded features
├── {subject}.sqlite × 168 subjects               # SQLite databases
└── patient_labels_updatedSept302025.xlsx         # Disease annotations

INPUT DIMENSIONS:
├── Subjects: 168 (across 5 diagnostic categories)
├── Initial features: 326 mitochondrial CellProfiler features
└── Cells per subject: Variable (after QC filtering)
```

### Processing Pipeline

#### Step 1: Data Loading & Merging (Lines 83-112)

```python
OPERATIONS:
1. pd.read_pickle() → df_1_0_0, df_1_0
2. read_single_cell_sql.readSingleCellData_sqlalch() → df_new
3. pd.merge() with disease_labels → df_1
4. Filter to common_subjects

OUTPUT: df_1 (merged single-cell data with labels)
```

#### Step 2: Quality Control (Lines 114-168)

```python
QC FILTERS APPLIED:
1. extract_cpfeature_names()
   → cp_features_analysis_0 (326 features)

2. handle_nans(thrsh_null_ratio=0.05, thrsh_std=0.001)
   → Drop features with >5% nulls or std < 0.001
   → cp_features_analysis (cleaned)

3. Border cell removal (borderLength=200px)
   Filter: !(Nuclei_Location_Center_X > 1188 OR < 200 OR
            Nuclei_Location_Center_Y > 840 OR < 200)
   → df_1_centCells

4. Segmentation quality filters:
   - Cells2Nuclei_MajorAxisLengthRatio > 2
   - Cells2Nuclei_AreaShapeRatio > 5

5. Intensity artifact removal:
   - Cells_Intensity_MeanIntensity_Actin < 0.5
   - Nuclei_Intensity_MeanIntensity_Actin < 0.55

OUTPUT: df_1 (clean single-cell data)
DIMENSIONS: Reduced from ~millions to cleaned cell count
```

#### Step 3: Feature Engineering (Lines 181-195)

```python
TRANSFORMATIONS:
1. StandardScaler().fit_transform()
   → df_1_scaled (z-scored features per feature)

2. find_correlation(threshold=0.7, remove_negative=True)
   → Remove highly correlated features
   → cp_features_collcorr

3. Label standardization:
   "MDD or Dep" → "MDD"

OUTPUT: df_1_scaled with standardized features and labels
```

#### Step 4: Representative Feature Selection (Lines 215-266)

```python
HIERARCHICAL CLUSTERING:
1. Filter to mito_features:
   - Contains "mito_" OR "RadialDistribution"
   - NOT "Actin" NOT "Correlation"
   → mito_features (subset of cp_features_analysis)

2. Compute correlation matrix:
   distance_matrix = 1 - |correlation|

3. hierarchy.ward(squareform(distance_matrix))
   → dist_linkage

4. hierarchy.fcluster(dist_linkage, y_ax_thsh=1.1)
   → cluster_ids

5. Select medoid from each cluster using hdmedians.medoid()
   → cluster_rep_feats (45 representative features)

6. Variance Inflation Factor (VIF) check:
   → All VIF < 10 (multicollinearity controlled)

OUTPUT: cluster_rep_feats (45 features)
SAVED: results/dendrogram_mito.pdf
```

#### Step 5: Machine Learning Classification (Lines 315-464)

```python
FOR EACH disease_label IN [psychosis, BP, SZ, SZA, MDD]:

    PREPARE DATA:
    1. If label == "psychosis":
          Combine BP, SZ, SZA → "psychosis"
    2. Filter to [Control, label]
    3. groupby([subject, label]).quantile(0.5)  # Per-patient median

    CROSS-VALIDATION:
    StratifiedKFold(n_splits=10)

    FOR EACH fold:
        1. LogisticRegression(solver='liblinear', class_weight='balanced')
        2. clf.fit(X_train, y_train)
        3. Calculate metrics:
           - balanced_accuracy_score()
           - roc_auc_score()

        4. Null distribution (n_permutations=10):
           - Shuffle y_train
           - Fit and evaluate
           → null_fold_aucs

        5. Statistical test:
           ttest_ind(auc_scores, combined_null_aucs)
           → p_value

    6. Feature importance:
       - Model coefficients: clf.coef_[0]
       - permutation_importance(n_repeats=10)

OUTPUT:
├── table_pval_ls (AUC results per fold)
├── mi_df_all (mutual information scores)
└── results/Figure3.pdf (ROC-AUC comparison)

KEY RESULTS:
- Psychosis: AUC = 0.68 (p < 0.05)
- BP: AUC = 0.72 (p < 0.05)
- SZ: AUC = 0.55 (marginal)
- SZA: AUC = 0.65 (p < 0.05)
- MDD: AUC = 0.55 (not significant)
```

#### Step 6: MITO-SLOPE Calculation (Lines 854-1050)

```python
ALGORITHM (find_end_slope2):

1. Define target_columns:
   "Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16"
   for i in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
   → 12 radial bins from middle of cell to edge

2. Per-subject aggregation:
   df_1_avg_persub = groupby([subject, label])[target_columns].mean()

3. Control baseline:
   ctrl = df_1_avg_persub[label == 'Control'][target_columns].mean()

4. Control-subtracted pattern:
   diff_pattern = subject_pattern - ctrl

5. Apply smoothing:
   savgol_filter(data, window_length=5, polyorder=3)

6. Find last extremum:
   min_max_indc = [argmax(data), argmin(data)]
   Filter to indices < len(data) - 2  # Exclude edge artifacts
   last_peak_ind = max(valid_indices)

7. Calculate slope:
   slope = (data[B-1] - data[last_peak_ind]) / (B-1 - last_peak_ind)

   Units: z-score change per bin

   Interpretation:
   - slope > 0: Peripheral enrichment (mitochondria at edges)
   - slope < 0: Central shift (mitochondria near nucleus)

STATISTICAL TESTING:
FOR EACH disease_label:
    ttest_ind(Control_slope, disease_slope, equal_var=False)
    → p_value, t_statistic

RESULTS:
- Psychosis: slope = -0.XX, p = 0.0069 ✓
- BP: slope = -0.XX, p = 0.0061 ✓
- SZ: slope = -0.XX, p = 0.017 ✓
- SZA: slope = -0.XX, p = 0.56 (trend)
- MDD: slope = +0.XX, p = 0.18 (opposite trend)

OUTPUT:
├── df_1_avg_persub with [last_peak_loc, slope] columns
├── results/Figure4b.pdf (Radial distribution patterns)
└── results/Figure4c.pdf (MITO-SLOPE boxplots)
```

### Output Files

```
GENERATED FILES:
results/
├── Figure3.pdf              # ML classification ROC-AUC boxplots
├── SuppFigure2.pdf          # CP-LR vs CNN comparison
├── Figure4b.pdf             # Radial distribution patterns by diagnosis
├── Figure4c.pdf             # MITO-SLOPE boxplots by diagnosis
└── dendrogram_mito.pdf      # Feature clustering dendrogram

DATA PRODUCTS:
└── (No intermediate data files saved - figures only)
```

---

## Notebook 2.0: Virtual Screening Pipeline

**File:** `notebooks/2.0-mh-virtual-screen.py`

This notebook implements a 4-step pipeline to screen public datasets for perturbations affecting mitochondrial localization.

### Configuration: Dataset-Specific Parameters (Lines 120-206)

```python
DATASET CONFIGURATIONS (ds_info_dict):

datasets = {
    'lincs': {
        'profiles_path': '...DrugRepurposing.../backend/',
        'meta_cols': ['Metadata_broad_sample', 'Metadata_dose_recode',
                      'Metadata_pert_id', 'Metadata_InChIKey14',
                      'Metadata_pert_type', 'Metadata_moa', 'Metadata_target',
                      'Metadata_pert_id_dose', 'Metadata_pert_name'],
        'pert_col': 'Metadata_pert_id_dose',
        'target_features_list': ['slope']
    },
    'CDRP': {
        'profiles_path': '...bioactivecompoundprofiling/backend/',
        'meta_cols': ['Metadata_broad_sample', 'Metadata_mmoles_per_liter2',
                      'Metadata_pert_id', 'Metadata_Sample_Dose', 'Metadata_moa'],
        'pert_col': 'Metadata_Sample_Dose',
        'target_features_list': ['slope']
    },
    'jump_orf': {
        'profiles_path': '...cpg0016-jump/source_4/backend/',
        'meta_cols': ['Metadata_Symbol', 'Metadata_broad_sample', 'Metadata_JCP2022'],
        'pert_col': 'Metadata_JCP2022',
        'target_features_list': ['slope']
    },
    'jump_crispr': {
        'profiles_path': '...cpg0016-jump/source_13/backend/',
        'meta_cols': ['Metadata_NCBI_Gene_ID', 'Metadata_Symbol', 'Metadata_JCP2022'],
        'pert_col': 'Metadata_JCP2022',
        'target_features_list': ['slope']
    },
    'jump_compound': {
        'profiles_path': '...cpg0016-jump/source_13/backend/',
        'meta_cols': ['Metadata_InChIKey', 'Metadata_InChI', 'Metadata_JCP2022'],
        'pert_col': 'Metadata_JCP2022',
        'target_features_list': ['slope']
    },
    'taorf': {
        'profiles_path': '...cpg0017-rohban-pathways/backend/',
        'meta_cols': ['Metadata_gene_name', 'Metadata_pert_name',
                      'Metadata_broad_sample', 'Metadata_moa'],
        'pert_col': 'Metadata_broad_sample',
        'target_features_list': ['slope']
    }
}
```

### Step 1: Metadata Preprocessing (Lines 115-234)

```python
PURPOSE: Standardize metadata format across datasets

FOR EACH dataset:

    INPUT: Raw metadata CSV

    OPERATIONS:
    1. Read dataset-specific metadata files
    2. Create unified columns:
       - batch_plate = Batch + "-" + Metadata_Plate
       - ctrl_well = boolean (dataset-specific control identifiers)
    3. Add Metadata_pert_type if not present

    CONTROL IDENTIFIERS BY DATASET:
    - LINCS: Metadata_pert_type == "control"
    - CDRP: Metadata_pert_type == "control"
    - JUMP-ORF: Metadata_Symbol in ["LacZ", "BFP", "HcRed", "LUCIFERASE"]
    - JUMP-CRISPR: Metadata_Symbol == "non-targeting"
    - JUMP-Compound: (compound-specific controls)
    - TA-ORF: Metadata_gene_name in ["LacZ", "Luciferase"]

    OUTPUT: workspace/metadata/preprocessed/annot_{dataset}.csv

    COLUMNS:
    ├── Metadata_Plate
    ├── Metadata_Well
    ├── Metadata_Batch
    ├── batch_plate
    ├── ctrl_well
    ├── Metadata_pert_type
    └── [dataset-specific metadata columns]

EXAMPLE OUTPUT STRUCTURE:
annot_jump_orf.csv:
    Metadata_Plate  Metadata_Well  Metadata_Batch  batch_plate  ctrl_well  Metadata_Symbol  Metadata_JCP2022
    BR00116991      A01            2021_04_26_...  2021...991   False      TP53             JCP2022_...
    BR00116991      A02            2021_04_26_...  2021...991   True       LacZ             JCP2022_...
```

### Step 2: Identify Orthogonal Features (Lines 762-916)

```python
PURPOSE: Find features uncorrelated with MITO-SLOPE for use as specificity filter

FOR EACH dataset:

    INPUT:
    ├── annot_{dataset}.csv (from Step 1)
    └── Backend profiles (well-level or sampled)

    PROCESS:

    1. Load profiles:
       IF dataset == 'lincs':
           # LINCS lacks well-level profiles
           sample_single_cells_from_sql(n_rand_ims=100 per plate)
       ELSE:
           read_per_well_data() OR read_per_well_data_csvs()

    2. Feature preprocessing:
       extract_cpfeature_names() → cp_features, cp_features_analysis_0
       handle_nans(thrsh_null_ratio=0.05, thrsh_std=0.001)
       → cp_features_analysis

    3. Merge with metadata:
       pd.merge(df_sag, annot, on=common_cols)

    4. Standardization:
       standardize_per_catX(df_sag, 'batch_plate', cp_features_analysis)
       → Per-plate z-scoring

    5. Calculate radial pattern slope:
       target_columns = ['Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16'
                        for i in range(5, 17)]  # Bins 5-16

       control_mean = df_sag[ctrl_well].groupby('batch_plate')[target_columns].mean()
       diff_pattern = df_sag[target_columns] - control_mean

       slope = np.apply_along_axis(find_end_slope2, 1, diff_pattern)
       df_sag['slope'] = slope[:, 1]

    6. Find orthogonal features:
       target_features_list = ['slope',
                              'Cytoplasm_RadialDistribution_MeanFrac_mito_tubeness_16of16',
                              'Cells_RadialDistribution_MeanFrac_mito_tubeness_16of16']

       uncorr_with_tfs = set()
       FOR tfeat IN target_features_list:
           corr_matrix = df_sag[cp_features_analysis].corrwith(df_sag[tfeat]).abs()
           uncorr_with_tfs &= set(features where |corr| < 0.1)

       # Further filter to exclude radial and correlation features
       uncorr_with_tfs = [f for f in uncorr_with_tfs
                         if "RadialDistribution" not in f
                         and "_Correlation" not in f]

    7. Remove redundancy within orthogonal set:
       similar_fs_2remove = find_correlation(df_sag[uncorr_with_tfs],
                                            threshold=0.6,
                                            remove_negative=True)
       uncorr_feats_condese = uncorr_with_tfs - similar_fs_2remove

    OUTPUT:
    └── results/target_pattern_orth_features_lists/{dataset}_2.csv
        Single column: orth_fs

    ORTHOGONAL FEATURE COUNTS:
    - LINCS: 7 features (manually curated due to sampling)
    - CDRP: 23 features
    - JUMP-ORF: 40 features
    - JUMP-CRISPR: 14 features
    - JUMP-Compound: Similar to CRISPR
    - TA-ORF: 14 features

EXAMPLE ORTHOGONAL FEATURES:
['Nuclei_AreaShape_FormFactor',
 'Nuclei_AreaShape_Eccentricity',
 'Cells_AreaShape_Solidity',
 'Cells_Intensity_MaxIntensity_Mito',
 'Cells_AreaShape_Eccentricity',
 'Cytoplasm_AreaShape_MaxFeretDiameter',
 'Nuclei_Texture_AngularSecondMoment_DNA_8_45']
```

### Step 3: Per-Site Aggregation (Lines 948-1003)

```python
PURPOSE: Create site-level profiles for target + orthogonal features

FOR EACH dataset:

    INPUT:
    ├── annot_{dataset}.csv
    ├── Orthogonal features list:
    │   - fibroblast_derived.csv (preferred) OR
    │   - {dataset}.csv (dataset-specific)
    └── Backend SQLite databases

    FEATURE LIST CONSTRUCTION:
    IF dataset == 'lincs':
        # Use manually curated set
        uncorr_feats_condese = ['Nuclei_AreaShape_FormFactor', ...]
    ELSE:
        # Use fibroblast-derived for consistency
        uncorr_feats_condese = read_csv('fibroblast_derived.csv')['orth_fs']

    radial_meanFrac_features = [
        'Cells_RadialDistribution_MeanFrac_mito_tubeness_{i}of16'
        for i in range(1, 17)  # ALL 16 bins
    ]

    feature_list2 = uncorr_feats_condese + radial_meanFrac_features

    PROCESS (form_per_site_aggregated_profiles):

    FOR EACH batch IN annot['Batch'].unique():

        FOR EACH plate IN batch:

            1. Read single-cell data:
               fileName = profiles_path + batch + "/" + plate + "/" + plate + ".sqlite"
               sc_df = readSingleCellData_sqlalch_features_subset(
                   fileName, feature_list2
               )

            2. Aggregate to site level:
               per_site_aggregate = sc_df.groupby(['Metadata_Well', 'Metadata_Site'])
                                         .mean(numeric_only=True)[feature_list2 + ['Count_Cells']]
                                         .reset_index()

            3. Add metadata:
               per_site_aggregate['Metadata_Batch'] = batch
               per_site_aggregate['Metadata_Plate'] = plate

        # Concatenate all plates in batch
        df_sag = pd.concat(plate_list, axis=0)

        # Save batch-level aggregation
        OUTPUT: per_site_aggregated_profiles_newpattern_2/
                {dataset}/{batch}_site_agg_profiles.csv.gz

DIMENSIONS:
- Rows: N_wells × N_sites_per_well (typically 4-9 sites/well)
- Columns: ~20-50 features + metadata
```

### Step 4: Statistical Testing (Lines 1021-1379)

```python
PURPOSE: Test each perturbation vs controls for MITO-SLOPE effect

INPUT:
├── annot_{dataset}.csv
├── Orthogonal features list (fibroblast_derived.csv)
└── All batch_site_agg_profiles.csv.gz files

LOADING & PREPROCESSING (Lines 1078-1169):

1. Load all per-site profiles:
   per_site_df_ls = []
   FOR EACH batch:
       fileNameToSave = per_site_profiles_path + dataset + "/" +
                       batch + "_site_agg_profiles.csv.gz"
       per_site_df_b = pd.read_csv(fileNameToSave)

       # Track low-variance features per plate
       cols2remove_lowVars_eachPlate += features where std < 0.001

       per_site_df_ls.append(per_site_df_b)

   per_site_df = pd.concat(per_site_df_ls, axis=0)

2. Clean features:
   handle_nans(per_site_df, target_columns + uncorr_feats_condese,
              thrsh_null_ratio=0.05, thrsh_std=0.001)

   uncorr_feats_cond = uncorr_feats_condese - cols2remove_lowVars_eachPlate

3. Merge with metadata:
   merge_how = 'inner' if dataset in ['jump_crispr', 'jump_compound'] else 'left'
   per_site_df = pd.merge(per_site_df, annot, how=merge_how, on=common_cols)

4. Standardization (per-plate z-scoring):
   per_site_df = standardize_per_catX(per_site_df, 'batch_plate',
                                      target_columns + uncorr_feats_cond)

5. Calculate MITO-SLOPE:
   # Get control mean pattern per plate
   control_df_perplate = per_site_df[ctrl_well].groupby('batch_plate')[target_columns].mean()

   # Subtract control mean within each plate
   df_rep_level_scaled_meanSub = per_site_df.groupby('batch_plate')[target_columns]
                                            .apply(subtract_control)

   # Apply find_end_slope2 to each row
   peak_slope = np.apply_along_axis(find_end_slope2, 1,
                                    df_rep_level_scaled_meanSub.values)

   per_site_df[['last_peak_loc', 'slope']] = peak_slope

   # Standardize slope
   per_site_df = standardize_per_catX(per_site_df, 'batch_plate',
                                      ... + ['last_peak_loc', 'slope'])

STATISTICAL TESTING (Lines 1191-1379):

pert_col = ds_info_dict[dataset]['pert_col']
meta_cols = ds_info_dict[dataset]['meta_cols']

results = annot[meta_cols].drop_duplicates().reset_index()

# Get perturbations (exclude controls)
IF 'Metadata_pert_type' in columns:
    perts = per_site_df[Metadata_pert_type in ['trt', 'Treated']][pert_col].unique()
ELSE:
    perts = per_site_df[~ctrl_well][pert_col].unique()

FOR EACH pert IN perts:

    # Filter to this perturbation
    per_site_df_pert = per_site_df[pert_col == pert]

    # Get plates where this perturbation exists (with >1 site)
    plates_pert = per_site_df_pert.groupby('batch_plate')
                                  .filter(lambda x: len(x) > 1)['batch_plate']
                                  .unique()

    # Initialize storage for per-plate statistics
    pert_cell_count_perSite_all_plates = []
    pert_pvals_all_plates = np.full((len(plates_pert), 6), np.nan)
    pert_tvals_all_plates = np.full((len(plates_pert), 4), np.nan)
    peak_slope_all_plates = np.full((len(plates_pert), 2), np.nan)

    FOR pi, plate IN enumerate(plates_pert):

        # Filter to this plate
        per_site_df_pert_plate = per_site_df_pert[batch_plate == plate]
        control_df = per_site_df[ctrl_well & batch_plate == plate]

        # Cell count (for filtering later)
        pert_cell_count_perSite_all_plates[pi] = per_site_df_pert_plate['Count_Cells'].mean()

        # TEST 1: Two-sample t-test on slope
        test_res = ttest_ind(per_site_df_pert_plate['slope'],
                            control_df['slope'],
                            equal_var=False)
        pert_pvals_all_plates[pi, 2] = test_res.pvalue
        pert_tvals_all_plates[pi, 2] = test_res.statistic

        # Calculate Cohen's d (effect size)
        cohend = cohens_d(per_site_df_pert_plate['slope'],
                         control_df['slope'])
        pert_tvals_all_plates[pi, 3] = cohend  # MAIN RANKING METRIC

        # Calculate standardized p-value
        degfree = len(pert_slope) + len(ctrl_slope) - 2
        z_score = t_to_z(test_res.statistic, degfree)
        std_p_val = z_to_p(z_score)
        pert_pvals_all_plates[pi, 3] = std_p_val

        # TEST 2: Hotelling's T² on full radial pattern
        statistic, p_value, p_value_std = TwoSampleT2Test(
            control_df[target_columns],
            per_site_df_pert_plate[target_columns]
        )
        pert_pvals_all_plates[pi, 0] = p_value
        pert_tvals_all_plates[pi, 0] = statistic
        pert_pvals_all_plates[pi, 4] = p_value_std

        # TEST 3: Hotelling's T² on orthogonal features
        statistic, p_value, p_value_std = TwoSampleT2Test(
            control_df[uncorr_feats_cond],
            per_site_df_pert_plate[uncorr_feats_cond]
        )
        pert_pvals_all_plates[pi, 1] = p_value
        pert_tvals_all_plates[pi, 1] = statistic
        pert_pvals_all_plates[pi, 5] = p_value_std

        # Store median slope values
        peak_slope_all_plates[pi, :] = per_site_df_pert_plate[['last_peak_loc', 'slope']].median()

    # SELECT REPRESENTATIVE PLATE (median by Cohen's d)
    med_t = np.nanpercentile(pert_tvals_all_plates, 50, axis=0, interpolation='nearest')
    median_selection_ind = 3  # Use d_slope (Cohen's d) as selection criterion

    IF ~np.isnan(med_t[median_selection_ind]):
        median_t_indx_val = argwhere(pert_tvals_all_plates[:, 3] == med_t[3])[0][0]
        median_t_indx = [median_t_indx_val] * 4
        median_t_indx_p = [median_t_indx_val] * 6
    ELSE:
        median_t_indx = [np.nan] * 4
        median_t_indx_p = [np.nan] * 6

    # Store results from median plate
    results.loc[pert_col == pert, 'Count_Cells_avg'] = mean(cell_counts)
    results.loc[pert_col == pert,
               ['p_target_pattern', 'p_orth', 'p_slope',
                'p_slope_std', 'p_pattern_std', 'p_orth_std']] = [
        pert_pvals_all_plates[median_t_indx_p[i], i] for i in range(6)
    ]
    results.loc[pert_col == pert,
               ['t_target_pattern', 't_orth', 't_slope', 'd_slope']] = [
        pert_tvals_all_plates[median_t_indx[i], i] for i in range(4)
    ]
    results.loc[pert_col == pert,
               ['last_peak_ind', 'slope']] = np.nanmedian(peak_slope_all_plates, axis=0)

# Sort by Cohen's d (effect size)
results = results.sort_values(by='slope', ascending=False)

OUTPUT: results/virtual_screen/{dataset}_results_pattern_aug_070624.csv

COLUMNS:
├── pert_col (dataset-specific: Metadata_JCP2022, Metadata_pert_id_dose, etc.)
├── meta_cols (dataset-specific metadata)
├── Count_Cells_avg (mean cell count across plates)
├── p_target_pattern (p-value for Hotelling's T² on full pattern)
├── p_orth (p-value for Hotelling's T² on orthogonal features)
├── p_slope (p-value for t-test on slope)
├── p_slope_std (standardized p-value for slope)
├── p_pattern_std (standardized p-value for pattern)
├── p_orth_std (standardized p-value for orthogonal)
├── t_target_pattern (t-statistic for pattern)
├── t_orth (t-statistic for orthogonal)
├── t_slope (t-statistic for slope)
├── d_slope (Cohen's d for slope) ← PRIMARY RANKING METRIC
├── last_peak_ind (location of last extremum in pattern)
└── slope (median MITO-SLOPE value)
```

### Output Files

```
GENERATED FILES:
workspace/
├── metadata/preprocessed/
│   ├── annot_lincs.csv
│   ├── annot_CDRP.csv
│   ├── annot_jump_orf.csv
│   ├── annot_jump_crispr.csv
│   ├── annot_jump_compound.csv
│   └── annot_taorf.csv
│
├── results/target_pattern_orth_features_lists/
│   ├── fibroblast_derived.csv       [Preferred - used across datasets]
│   ├── lincs_2.csv                  [7 features]
│   ├── CDRP_2.csv                   [23 features]
│   ├── jump_orf_2.csv               [40 features]
│   ├── jump_crispr_2.csv            [14 features]
│   └── taorf_2.csv                  [14 features]
│
├── per_site_aggregated_profiles_newpattern_2/
│   ├── lincs/{batch}_site_agg_profiles.csv.gz
│   ├── CDRP/{batch}_site_agg_profiles.csv.gz
│   ├── jump_orf/{batch}_site_agg_profiles.csv.gz
│   ├── jump_crispr/{batch}_site_agg_profiles.csv.gz
│   ├── jump_compound/{batch}_site_agg_profiles.csv.gz
│   └── taorf/{batch}_site_agg_profiles.csv.gz
│
└── results/virtual_screen/
    ├── lincs_results_pattern_aug_070624.csv
    ├── CDRP_results_pattern_aug_070624.csv
    ├── jump_orf_results_pattern_aug_070624.csv
    ├── jump_crispr_results_pattern_aug_070624.csv
    ├── jump_compound_results_pattern_aug_070624.csv
    └── taorf_results_pattern_aug_070624.csv

DIMENSIONS:
- LINCS: ~9,426 rows (1,571 compounds × 6 doses)
- CDRP: ~30,430 rows (compounds × doses)
- JUMP-ORF: ~12,609 rows (genes)
- JUMP-CRISPR: ~7,975 rows (genes)
- JUMP-Compound: ~116,000 rows (compounds)
- TA-ORF: ~323 rows (genes)
```

---

## Notebook 2.1: Gene Set Enrichment Analysis

**File:** `notebooks/2.1-mh-set-enrichment-analysis.py`

### Input Data

```python
INPUT FILES:
└── results/virtual_screen/{dataset}_results_pattern_aug_070624.csv
    (One file per dataset from 2.0 Step 4)

EXTERNAL GENE SET LIBRARIES (via blitzgsea):
├── GO_Biological_Process_2023
├── GO_Molecular_Function_2023
├── GO_Cellular_Component_2023
├── KEGG_2021_Human
├── WikiPathway_2021_Human
├── PFOCR_Pathways_2023
├── OMIM_Disease
├── OMIM_Expanded
├── Human_Phenotype_Ontology
├── MGI_Mammalian_Phenotype_Level_4_2021
├── Proteomics_Drug_Atlas_2023
└── GWAS_Catalog_2023

EXTERNAL DRUG ANNOTATIONS:
└── workspace/metadata/CompoundClusters08202104.xlsx
    (Rakesh's drug mechanism-of-action categories)

CUSTOM GENE LISTS (from Cuperfain et al. 2018 Table 3):
├── Fusion genes: [AFG3L2, BAK1, BAX, BNIP3, CHCHD3, DNM1L, ...]
├── Fission genes: [BNIP3, COX10, DHODH, DNM1L, DNM2, DNM3, ...]
├── Transport genes: [CLUH, DNM1L, KIF1B, LRPPRC, LRRK2, ...]
└── Mitophagy genes: [ATPIF1, BNIP3, BNIP3L, CISD2, DNM1L, ...]
```

### Processing Pipeline

#### Configuration (Lines 92-109)

```python
# BH-corrected critical values for orthogonal filter (from 2.2)
orth_bh_corrected_critical_dict = {
    'taorf': 0.0329,
    'lincs': 0.0483,
    'CDRP': 0.0481,
    'jump_orf': 0.0375,
    'jump_crispr': 0.0479,
    'jump_compound': 0.0475
}

# Dataset-specific column mappings
datasets_info_dict = {
    'taorf': {
        'key_col': 'Metadata_gene_name',
        'key_sample_id_col': 'Metadata_broad_sample'
    },
    'jump_orf': {
        'key_col': 'Metadata_Symbol',
        'key_sample_id_col': 'Metadata_broad_sample'  # OR Metadata_JCP2022
    },
    'jump_crispr': {
        'key_col': 'Metadata_Symbol',
        'key_sample_id_col': 'Metadata_JCP2022'
    },
    'jump_compound': {
        'key_col': 'Metadata_InChIKey',
        'key_col_ref_set': 'Metadata_JCP2022'
    },
    'lincs': {
        'key_col': 'Metadata_pert_name',
        'key_col_ref_set': 'Metadata_pert_id_dose'
    },
    'CDRP': {
        'key_col': 'Metadata_pert_id',
        'key_col_ref_set': 'Metadata_Sample_Dose'
    }
}
```

#### Analysis 1: GO/KEGG/WikiPathways Enrichment (Lines 128-210)

```python
FOR dataset IN ['taorf', 'jump_orf', 'jump_crispr']:

    1. LOAD & FILTER:
       res_df = read_csv(virtual_screen/{dataset}_results_pattern_aug_070624.csv)
       res_df = res_df[~isnull(key_col) & ~isnull(d_slope)]
       res_df = groupby(key_col).median()  # Per-gene median

       IF cell_count_filter_enabled:
           res_df = res_df[Count_Cells_avg > quantile(0.1)]

       IF orth_filter_enabled:
           res_df = res_df[p_orth_std > orth_bh_corrected_critical]

    2. PREPARE SIGNATURE:
       sig_df = res_df[[key_col, d_slope]]
               .rename(columns={key_col: 0, d_slope: 1})
               .sort_values(by=1)

       # Signature format for blitzGSEA:
       # Column 0: Gene identifier
       # Column 1: Ranking metric (d_slope)

    3. RUN ENRICHMENT:
       FOR database_str IN database_ls:

           library = blitz.enrichr.get_library(database_str)
           result = blitz.gsea(sig_df, library, min_size=4, seed=1)

           top_res_df = result[result['fdr'] < 0.05]

           IF len(top_res_df) > 0:
               # Generate visualization
               fig_table = top_table_custom(sig_df, library, top_res_df,
                                           set_label="Gene Set",
                                           n=top_res_df.shape[0])
               plt.show()

               PRINT: "Dataset: {dataset}, {top_res_df.shape[0]} / {result.shape[0]} categories enriched"

OUTPUT: Console output + figures (not saved to files)
```

#### Analysis 2: Mitochondrial Dynamics Genes (Lines 222-487)

```python
PURPOSE: Test enrichment of custom mitochondrial gene sets from literature

GENE LISTS (from Cuperfain et al. 2018, Table 3):
Fusion_ls = ['AFG3L2', 'BAK1', 'BAX', 'BNIP3', 'CHCHD3', 'DNM1L',
             'FIS1', 'GDAP1', 'MFF', 'MFN1', 'MFN2', 'MIEF1', 'MIEF2',
             'MUL1', 'OMA1', 'OPA1', 'PLD6', 'PRKN', 'STOML2', 'USP30']

Fission_ls = ['BNIP3', 'COX10', 'DHODH', 'DNM1L', 'DNM2', 'DNM3',
              'FIS1', 'GDAP1', 'LRRK2', 'MARCH5', 'MFF', 'MIEF1',
              'MIEF2', 'MTFP1', 'MTFR1', 'MTFR1L', 'MTFR2', 'MUL1',
              'MYO19', 'OPA1', 'PINK1', 'PRKN', 'SLC25A46']

Transport_ls = ['CLUH', 'DNM1L', 'KIF1B', 'LRPPRC', 'LRRK2', 'MAIP1',
                'MAP1S', 'MFN1', 'MFN2', 'MGARP', 'MST01', 'MUL1',
                'OPA1', 'RHOT1', 'RHOT2', 'SYNJ2BP', 'TRAK1', 'TRAK2', 'UBB']

Mitophagy_ls = ['ATPIF1', 'BNIP3', 'BNIP3L', 'CISD2', 'DNM1L', 'FIS1',
                'FUNDC1', 'FUNDC2', 'HK2', 'HTRA2', 'MFN2', 'MUL1',
                'PARK7', 'PINK1', 'SQSTM1', 'TOMM7', 'TSPO', 'VDAC1']

table3_list = list(set(Fusion_ls + Fission_ls + Transport_ls + Mitophagy_ls))

FOR dataset IN ['jump_orf', 'jump_crispr']:

    key_col = datasets_info_dict[dataset]['key_col']
    key_sample_id_col = datasets_info_dict[dataset]['key_sample_id_col']

    # Load and aggregate by sample ID (not gene symbol)
    res_df = read_csv(virtual_screen/{dataset}_results_pattern_aug_070624.csv)

    agg_dict = {col: 'median' if numeric else 'first'
                for col in columns if col != key_sample_id_col}

    res_df = groupby(key_sample_id_col).agg(agg_dict).reset_index()

    # Map gene lists to sample IDs
    Transport_ls_brd = res_df.loc[res_df[key_col].isin(Transport_ls),
                                  key_sample_id_col].unique().tolist()

    table3_ls_brd = res_df.loc[res_df[key_col].isin(table3_list),
                               key_sample_id_col].unique().tolist()

    # Create library with sample IDs (not gene symbols)
    library_table3 = {
        'Mitochondrial dynamics': table3_ls_brd,
        'Mitochondrial transport': Transport_ls_brd
    }

    # Prepare signature with sample IDs
    sig_df = res_df[[key_sample_id_col, d_slope]]
            .rename(columns={key_sample_id_col: 0, d_slope: 1})
            .sort_values(by=1)

    # Run enrichment
    result_table3set = blitz.gsea(sig_df, library_table3, min_size=1, seed=1)

    PRINT: result_table3set[['es', 'nes', 'pval', 'sidak', 'fdr', 'geneset_size']]

    fig_table = top_table_custom(sig_df, library_table3, result_table3set,
                                 set_label="Gene Set", n=2)
    plt.show()

EXPECTED RESULTS:
- JUMP-CRISPR: Mitochondrial dynamics genes should be enriched (NES < 0)
- JUMP-ORF: May show weaker or opposite enrichment
```

#### Analysis 3: Schizophrenia-Associated Genes (Lines 505-603)

```python
PURPOSE: Test genes from "Rare coding variants in ten genes confer substantial
         risk for schizophrenia" (Daly et al., Nature)

GENE LISTS:
SZ_genes_Daly = ['CACNA1G', 'GRIN2A', 'GRIA3', 'TRIO', 'SP4',
                 'RB1CC1', 'SETD1A', 'XPO7', 'CUL1', 'HERC1']

DD_genes_Daly = ['STAG1', 'ASH1L', 'ZMYM2', 'KDM6B', 'SRRM2',
                 'HIST1H1E', 'PMEPA1']

misc_genes = ['CACNA1A', 'AKAP11', 'SETD1A']

FOR dataset IN ['jump_orf', 'jump_crispr']:

    # Same aggregation approach as Analysis 2

    # Create library
    library_Daly = {}
    library_Daly_symbols = {}

    FOR gls_name, gls IN zip(['SZ genes', 'DD/ID genes', 'Misc genes'],
                             [SZ_genes_Daly, DD_genes_Daly, misc_genes]):

        # Map to sample IDs present in dataset
        library_Daly[gls_name] = res_df.loc[res_df[key_col].isin(gls),
                                            key_sample_id_col].unique()

        # Track which symbols were found
        library_Daly_symbols[gls_name] = res_df.loc[res_df[key_col].isin(gls),
                                                    key_col].unique()

    PRINT: library_Daly_symbols

    # Run enrichment
    result_Daly = blitz.gsea(sig_df, library_Daly, min_size=1, seed=1)

    PRINT: result_Daly[['es', 'nes', 'pval', 'sidak', 'fdr', 'geneset_size']]

    fig_table = top_table_custom(sig_df, library_Daly, result_Daly,
                                 set_label="Gene Set", n=2)
    plt.show()
```

#### Analysis 4: Drug MOA Categories - LINCS (Lines 608-803)

```python
PURPOSE: Test enrichment of Rakesh's manually curated drug categories in LINCS

INPUT:
└── workspace/metadata/CompoundClusters08202104.xlsx

DRUG CATEGORIES (examples):
- 5HT antagonist
- Acetylcholine receptor antagonist
- Adrenergic receptor agonist
- Adrenergic receptor antagonist
- Anesthetic - dissociative
- Antipsychotic
- Dopamine receptor agonist
- Dopamine receptor antagonist
- GABA receptor agonist
- Histamine receptor antagonist
- Lithium
- MAOI
- NMDA antagonist
- Norepinephrine reuptake inhibitor
- Opioid agonist
- Serotonin norepinephrine reuptake inhibitor
- Serotonin reuptake inhibitor
- Tricyclic Antidepressant
- And more...

PROCESS:

1. Load screening results:
   res_df_lincs = read_csv(virtual_screen/lincs_results_pattern_aug_070624.csv)
   res_df_lincs = res_df_lincs[~isnull(d_slope) & ~isnull(Metadata_pert_name)]
   res_df_lincs['Metadata_pert_name_lowercase'] = lower(Metadata_pert_name)

2. Optional filtering:
   IF cell_count_filter_enabled:
       res_df_lincs = res_df_lincs[Count_Cells_avg > quantile(0.1)]

   IF orth_filter_enabled:
       res_df_lincs = res_df_lincs[p_orth_std > orth_bh_corrected_critical]

3. Dose selection strategy:
   # Keep only dose with maximum absolute phenotype strength per compound
   res_df_lincs['d_slope_abs'] = abs(d_slope)
   abs_max_indices = groupby('Metadata_pert_name_lowercase')['d_slope_abs'].idxmax()
   res_df_lincs = res_df_lincs.loc[abs_max_indices]

4. Build library:
   drug_list_rakesh = read_excel(CompoundClusters08202104.xlsx)
   rakesh_cats = drug_list_rakesh.columns  # Drug categories

   lincs_prt_names = set(res_df_lincs['Metadata_pert_name'].str.lower())

   lib_rakesh = {}
   FOR category IN rakesh_cats:
       overlap_drugs = set(drug_list_rakesh[category].str.lower()) & lincs_prt_names
       overlap_drugs = {x for x in overlap_drugs if x == x}  # Remove NaNs

       IF len(overlap_drugs) > 0:
           # Map to reference set (pert_id or pert_id_dose)
           pert_id_dose_ls = res_df_lincs.loc[
               res_df_lincs['Metadata_pert_name_lowercase'].isin(overlap_drugs),
               key_col_ref_set
           ].unique()

           IF len(pert_id_dose_ls) > 1:
               lib_rakesh[category] = list(pert_id_dose_ls)

       PRINT: f"{category}: {len(overlap_drugs)} out of {total} exist in LINCS"

5. Run enrichment:
   sig_df_lincs = res_df_lincs[[key_col_ref_set, d_slope]]
                  .rename(columns={key_col_ref_set: 0, d_slope: 1})
                  .sort_values(by=1)

   result_rakesh_enrich = blitz.gsea(sig_df_lincs, lib_rakesh,
                                     min_size=5, seed=1)

   fig_table = top_table_custom(sig_df_lincs, lib_rakesh, result_rakesh_enrich,
                                set_label="Drug Set", n=result.shape[0])
   plt.show()

EXPECTED RESULTS:
- Most psychiatric drug classes should NOT be enriched (novel mechanism)
- May see enrichment in related pathways (e.g., monoamine-related)
```

#### Analysis 5: MOA Enrichment for All Compounds (Lines 806-944)

```python
PURPOSE: Test enrichment of mechanism-of-action categories in compound datasets

FOR dataset IN ['lincs', 'CDRP']:

    1. LOAD & PROCESS:
       res_df = read_csv(virtual_screen/{dataset}_results_pattern_aug_070624.csv)
       res_df['Metadata_moa'] = str.lower(Metadata_moa)

       IF cell_count_filter_enabled:
           res_df = res_df[Count_Cells_avg > quantile(0.1)]

       IF orth_filter_enabled:
           res_df = res_df[p_orth_std > orth_bh_corrected_critical]

       res_df['d_slope_abs'] = abs(d_slope)

    2. OPTIONAL DOSE SELECTION:
       # Keep only max dose per compound
       abs_max_indices = groupby(key_col + '_lowercase')['d_slope_abs'].idxmax()
       res_df = res_df.loc[abs_max_indices]

    3. HANDLE MISSING MOAs:
       res_df['Metadata_moa'] = fillna('nan')

    4. BUILD MOA LIBRARY:
       # Extract unique MOA terms
       unq_moas0 = []
       all_moas = res_df['Metadata_moa'].unique()

       FOR moa IN all_moas:
           IF moa == moa:  # Not NaN
               unq_moas0 += moa.split('|')  # MOAs are pipe-separated

       unq_moas = list(set(unq_moas0))

       lib_moa = {}
       FOR moa_term IN unq_moas:
           pert_id_dose_ls = res_df.loc[
               res_df['Metadata_moa'].str.contains(moa_term),
               key_col_ref_set
           ].unique()

           IF len(pert_id_dose_ls) > 1:
               lib_moa[moa_term] = list(pert_id_dose_ls)

    5. RUN ENRICHMENT:
       sig_df = res_df[[key_col_ref_set, d_slope]]
               .rename(columns={key_col_ref_set: 0, d_slope: 1})
               .sort_values(by=1)

       min_n_compounds_per_cat = 10
       result_moa_enrich = blitz.gsea(sig_df, lib_moa,
                                      min_size=min_n_compounds_per_cat,
                                      seed=1)

       enriched_df = result_moa_enrich[fdr < 0.05]

       fig_table = top_table_custom(sig_df, lib_moa, result_moa_enrich,
                                    set_label="MOA", n=enriched_df.shape[0])
       plt.show()

       PRINT: f"{dataset} has {len(unq_moas)} MOA categories, "
              f"{result_moa_enrich.shape[0]} have min_size > {min_n_compounds_per_cat}"
       PRINT: f"{enriched_df.shape[0]} / {result_moa_enrich.shape[0]} enriched"
       PRINT: enriched_df.index
```

### Output

```
OUTPUT FORMAT: Console + Matplotlib figures (not saved to files)

TYPICAL OUTPUT:
┌────────────────────────────────────────────────────────────┐
│                      Gene Set Enrichment                   │
├────────────────────────────────────────────────────────────┤
│ FDR    NES    Gene Set                                     │
├────────────────────────────────────────────────────────────┤
│ 0.012  -2.45  Mitochondrial dynamics                       │
│ 0.034  -1.89  Mitochondrial transport                      │
└────────────────────────────────────────────────────────────┘

* jump_crispr: Out of 156 GO_Biological_Process_2023 categories,
               23 are enriched.
* jump_crispr: Out of 89 KEGG_2021_Human categories,
               12 are enriched.
```

---

## Notebook 2.2: Hit Filtering & Final Tables

**File:** `notebooks/2.2-mh-check-vs-lists.py`

### Input Data

```python
INPUT:
└── results/virtual_screen/{dataset}_results_pattern_aug_070624.csv
    × 6 datasets (from 2.0 Step 4)

CONSTANTS (from prior analysis):
└── results/target_pattern_orth_features_lists/fibroblast_derived.csv
```

### Processing Pipeline

#### Benjamini-Hochberg Correction Function (Lines 266-278)

```python
def bh_adjusted_critical_value(pvalues, fdr=0.05):
    """
    Calculate BH-corrected critical value for multiple hypothesis testing

    Returns the maximum p-value that passes BH correction at given FDR
    """
    sorted_pvalues = np.sort(pvalues)
    m = len(pvalues)
    ranks = np.arange(1, m + 1)
    critical_values = (ranks / m) * fdr

    below_threshold = sorted_pvalues <= critical_values

    IF np.any(below_threshold):
        adjusted_critical = sorted_pvalues[below_threshold].max()
    ELSE:
        adjusted_critical = np.nan

    RETURN adjusted_critical
```

#### Main Filtering Loop (Lines 285-442)

```python
# Configuration
sort_by_col = 'd_slope'                    # Cohen's d (effect size)
p_val_target_col = 'p_slope_std'           # Standardized p-value for slope
p_val_orth_col = 'p_orth_std'              # Standardized p-value for orthogonal

# Tracking dataframes
orth_bh_corrected_critical_dict = {}
target_bh_corrected_critical_dict = {}
filtering_stats = pd.DataFrame(columns=['raw', 'target feature significance',
                                       'orth filter', 'both filters'])

FOR dataset IN ['taorf', 'lincs', 'CDRP', 'jump_orf', 'jump_crispr', 'jump_compound']:

    # Get dataset configuration
    meta_cols = ds_info_dict[dataset]['meta_cols']
    pert_col = ds_info_dict[dataset]['pert_col']

    # LOAD RESULTS
    res_df = read_csv(virtual_screen/{dataset}_results_pattern_aug_070624.csv)
    res_df = res_df[~isnull(sort_by_col)]

    filtering_stats.loc[dataset, 'raw'] = res_df.shape[0]

    # OPTIONAL: Cell count filter (default: OFF in this notebook)
    IF cell_count_filter_enabled:
        res_df_hcc = res_df[Count_Cells_avg > quantile(0.1)]
    ELSE:
        res_df_hcc = res_df.copy()

    # CALCULATE BH-CORRECTED CRITICAL VALUES
    target_critical = bh_adjusted_critical_value(res_df_hcc[p_val_target_col],
                                                 fdr=0.05)
    orth_critical = bh_adjusted_critical_value(res_df_hcc[p_val_orth_col],
                                               fdr=0.05)

    orth_bh_corrected_critical_dict[dataset] = np.round(orth_critical, 5)
    target_bh_corrected_critical_dict[dataset] = np.round(target_critical, 5)

    PRINT: f"{dataset} orth_bh_corrected_critical = {orth_critical}"

    # FILTER 1: Target feature significance
    res_df_target_sig = res_df_hcc[p_val_target_col < target_critical]
    filtering_stats.loc[dataset, 'target feature significance'] = res_df_target_sig.shape[0]

    # FILTER 2: Orthogonal features non-significant (cells look normal)
    res_df_orth_filt = res_df_hcc[p_val_orth_col > orth_critical]
    filtering_stats.loc[dataset, 'orth filter'] = res_df_orth_filt.shape[0]

    # BOTH FILTERS
    res_df_both_filt = res_df_target_sig[res_df_target_sig[p_val_orth_col] > orth_critical]
    filtering_stats.loc[dataset, 'both filters'] = res_df_both_filt.shape[0]

    # PREPARE OUTPUT COLUMNS
    list_of_cols_2save = [pert_col, sort_by_col, p_val_target_col, p_val_orth_col,
                          't_orth', 'Count_Cells_avg'] + meta_cols

    # SORT BY EFFECT SIZE
    res_df_all = res_df.sort_values(by=sort_by_col)[list_of_cols_2save]
    res_df_orth = res_df_orth_filt.sort_values(by=sort_by_col)[list_of_cols_2save]
    res_df_both = res_df_both_filt.sort_values(by=sort_by_col)[list_of_cols_2save]

    # Remove duplicate columns
    res_df_all = res_df_all.loc[:, ~res_df_all.columns.duplicated(keep='last')]
    res_df_orth = res_df_orth.loc[:, ~res_df_orth.columns.duplicated(keep='last')]
    res_df_both = res_df_both.loc[:, ~res_df_both.columns.duplicated(keep='last')]

    # SAVE TO EXCEL (3 sheets)
    saveAsNewSheetToExistingFile(
        virtual_screen_results_202407/{dataset}_screen_results.xlsx,
        [res_df_all, res_df_orth, res_df_both],
        [dataset, dataset + '_orthfilt', dataset + '_bothfilt'],
        keep_index_column=False
    )

    # Store for downstream use
    list_of_res_df1[dataset] = res_df_orth
    list_of_res_df2[dataset] = res_df_both

    PRINT: f"{dataset}: original shape: {res_df.shape}, filtered shape: {res_df_both.shape}"

# PRINT SUMMARY TABLE
PRINT: filtering_stats
```

#### BH-Corrected Critical Values (Output)

```python
CALCULATED VALUES:

orth_bh_corrected_critical_dict = {
    'taorf': 0.03287,
    'lincs': 0.04829,
    'CDRP': 0.04812,
    'jump_orf': 0.03752,
    'jump_crispr': 0.04785,
    'jump_compound': 0.04746
}

target_bh_corrected_critical_dict = {
    'taorf': 0.00409,
    'lincs': 0.00761,
    'CDRP': 0.00924,
    'jump_orf': 0.00198,
    'jump_crispr': 0.007,
    'jump_compound': 0.00365
}

INTERPRETATION:
- Lower target critical values = more stringent filtering
- JUMP-ORF has strictest target filter (0.00198)
- All orth filters similar (~0.03-0.05)
```

#### Example Filtering Statistics

```
TYPICAL FILTERING PROGRESSION:

Dataset       | Raw    | Target Sig | Orth Filter | Both Filters
--------------|--------|------------|-------------|-------------
lincs         | 9,426  | 1,247      | 3,892       | 487
CDRP          | 30,108 | 3,821      | 12,456      | 1,892
jump_orf      | 12,609 | 892        | 5,234       | 421
jump_crispr   | 7,975  | 634        | 3,118       | 289
jump_compound | 116,234| 8,912      | 48,723      | 3,456
taorf         | 323    | 54         | 178         | 28

FILTERING RATIO:
- Target filter: ~10-15% of raw
- Orth filter: ~30-40% of raw
- Both filters: ~5-10% of raw (FINAL HIT LIST)
```

#### Specific Gene Validation (Lines 474-478)

```python
# Example: Check specific psychiatric genes
df = list_of_res_df['jump_orf']
genes_of_interest = ['GSK3B', 'GSK3A', 'CACNA1C']

validation_df = df.loc[
    df['Metadata_gene_name'].isin(genes_of_interest),
    ['Metadata_gene_name', 'Metadata_pert_name',
     'p_slope', 'p_slope_std', 'd_slope', 'slope']
]

PRINT: validation_df

EXPECTED:
   Metadata_gene_name  p_slope_std  d_slope  slope
   CACNA1C            0.0023       -0.45    -0.12
   GSK3A              0.0891        0.23     0.08
   GSK3B              0.0456        0.31     0.09
```

### Output Files

```
GENERATED FILES:
workspace/results/virtual_screen_results_202407/

├── lincs_screen_results.xlsx
│   ├── Sheet 1: lincs (all results, sorted by d_slope)
│   ├── Sheet 2: lincs_orthfilt (orthogonal filter only)
│   └── Sheet 3: lincs_bothfilt (both filters - FINAL HITS)
│
├── CDRP_screen_results.xlsx
│   ├── Sheet 1: CDRP
│   ├── Sheet 2: CDRP_orthfilt
│   └── Sheet 3: CDRP_bothfilt
│
├── jump_orf_screen_results.xlsx
│   ├── Sheet 1: jump_orf
│   ├── Sheet 2: jump_orf_orthfilt
│   └── Sheet 3: jump_orf_bothfilt
│
├── jump_crispr_screen_results.xlsx
│   ├── Sheet 1: jump_crispr
│   ├── Sheet 2: jump_crispr_orthfilt
│   └── Sheet 3: jump_crispr_bothfilt
│
├── jump_compound_screen_results.xlsx
│   ├── Sheet 1: jump_compound
│   ├── Sheet 2: jump_compound_orthfilt
│   └── Sheet 3: jump_compound_bothfilt
│
└── taorf_screen_results.xlsx
    ├── Sheet 1: taorf
    ├── Sheet 2: taorf_orthfilt
    └── Sheet 3: taorf_bothfilt

COLUMNS IN EXCEL FILES:
├── pert_col (dataset-specific perturbation identifier)
├── d_slope (Cohen's d - PRIMARY RANKING METRIC)
├── p_slope_std (standardized p-value for slope)
├── p_orth_std (standardized p-value for orthogonal features)
├── t_orth (t-statistic for orthogonal features)
├── Count_Cells_avg (mean cell count)
└── [dataset-specific metadata columns]

EXAMPLE HITS (Sheet 3 - both filters):

LINCS positive d_slope (increase mitochondria at edges):
   Metadata_pert_name      d_slope  p_slope_std  p_orth_std
   divalproex-sodium       0.87     0.0012       0.0654
   XL-147                  0.76     0.0023       0.0891
   tivozanib               0.68     0.0034       0.1123

LINCS negative d_slope (decrease mitochondria at edges):
   Metadata_pert_name      d_slope  p_slope_std  p_orth_std
   PF-02545920            -0.92     0.0008       0.0723
   citalopram             -0.81     0.0015       0.0834
   trapidil               -0.73     0.0021       0.0956
```

---

## Complete Data Dependency Graph

```
┌─────────────────────────────────────────────────────────────────┐
│                       EXTERNAL DATA SOURCES                      │
└───────┬─────────────────────────────────────────────────────────┘
        │
        ├─ Patient Fibroblasts (SQLite DBs, Pickles, Labels)
        └─ Public Cell Painting (Backend: SQLite/CSV/Parquet)
                                  Metadata: CSVs
                                  Gene sets: Enrichr libraries
        │
        ↓
┌───────────────────────────────────────────────────────────────────┐
│                         1.0 PATIENT ANALYSIS                      │
├───────────────────────────────────────────────────────────────────┤
│ Input:  SQLite DBs, Pickles, patient_labels.xlsx                 │
│ Process: QC → Feature Selection → ML → MITO-SLOPE                │
│ Output: Figures 3, 4 (PDFs)                                      │
│         [No intermediate data files]                              │
└───────────────────────────────────────────────────────────────────┘
        │
        │ (Defines MITO-SLOPE algorithm)
        │ (Identifies orthogonal features)
        ↓
┌───────────────────────────────────────────────────────────────────┐
│                   2.0 STEP 1: METADATA PREP                       │
├───────────────────────────────────────────────────────────────────┤
│ Input:  Raw metadata CSVs (dataset-specific)                     │
│ Output: workspace/metadata/preprocessed/annot_{dataset}.csv      │
└──────────────────────┬────────────────────────────────────────────┘
                       │
                       ↓
┌───────────────────────────────────────────────────────────────────┐
│                   2.0 STEP 2: FIND ORTHOGONAL                     │
├───────────────────────────────────────────────────────────────────┤
│ Input:  annot_{dataset}.csv                                      │
│         Backend profiles (well-level)                             │
│ Output: results/target_pattern_orth_features_lists/              │
│         ├── {dataset}_2.csv                                       │
│         └── fibroblast_derived.csv (consolidated)                │
└──────────────────────┬────────────────────────────────────────────┘
                       │
                       ↓
┌───────────────────────────────────────────────────────────────────┐
│                   2.0 STEP 3: PER-SITE AGG                        │
├───────────────────────────────────────────────────────────────────┤
│ Input:  annot_{dataset}.csv                                      │
│         fibroblast_derived.csv                                    │
│         Backend SQLite databases                                  │
│ Output: per_site_aggregated_profiles_newpattern_2/               │
│         {dataset}/{batch}_site_agg_profiles.csv.gz               │
└──────────────────────┬────────────────────────────────────────────┘
                       │
                       ↓
┌───────────────────────────────────────────────────────────────────┐
│                   2.0 STEP 4: STATISTICAL TEST                    │
├───────────────────────────────────────────────────────────────────┤
│ Input:  annot_{dataset}.csv                                      │
│         fibroblast_derived.csv                                    │
│         All {batch}_site_agg_profiles.csv.gz                     │
│ Output: results/virtual_screen/                                  │
│         {dataset}_results_pattern_aug_070624.csv                 │
└──────────────────────┬────────────────────────────────────────────┘
                       │
                       ├──────────────┬────────────────┐
                       ↓              ↓                ↓
            ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
            │ 2.1 ENRICHMENT│  │ 2.2 FILTERING│  │ MANUSCRIPT   │
            ├──────────────┤  ├──────────────┤  ├──────────────┤
            │ Input:       │  │ Input:       │  │ Input:       │
            │ - Results CSV│  │ - Results CSV│  │ - Figures    │
            │ - Gene sets  │  │              │  │ - Tables     │
            │              │  │ Process:     │  │              │
            │ Process:     │  │ - BH correct │  │ Output:      │
            │ - blitzGSEA  │  │ - 2 filters  │  │ - Paper      │
            │              │  │              │  │              │
            │ Output:      │  │ Output:      │  │              │
            │ - Console    │  │ - XLSX files │  │              │
            │ - Figures    │  │   (3 sheets) │  │              │
            └──────────────┘  └──────────────┘  └──────────────┘
```

---

## Key Intermediate Files Summary

```
workspace/
│
├── metadata/
│   ├── preprocessed/
│   │   ├── annot_lincs.csv              [2.0 Step 1 → Step 2,3,4]
│   │   ├── annot_CDRP.csv
│   │   ├── annot_jump_orf.csv
│   │   ├── annot_jump_crispr.csv
│   │   ├── annot_jump_compound.csv
│   │   └── annot_taorf.csv
│   │
│   ├── JUMP/                            [External - for metadata prep]
│   │   ├── plate.csv.gz
│   │   ├── well.csv.gz
│   │   ├── compound.csv.gz
│   │   ├── orf.csv.gz
│   │   └── crispr.csv.gz
│   │
│   └── CompoundClusters08202104.xlsx    [External - for 2.1 MOA analysis]
│
├── results/
│   ├── target_pattern_orth_features_lists/
│   │   ├── fibroblast_derived.csv       [2.0 Step 2 → Step 3,4]
│   │   ├── lincs_2.csv                  [Dataset-specific alternatives]
│   │   ├── CDRP_2.csv
│   │   └── ... (other datasets)
│   │
│   ├── virtual_screen/
│   │   ├── lincs_results_pattern_aug_070624.csv      [2.0 Step 4 → 2.1, 2.2]
│   │   ├── CDRP_results_pattern_aug_070624.csv
│   │   ├── jump_orf_results_pattern_aug_070624.csv
│   │   ├── jump_crispr_results_pattern_aug_070624.csv
│   │   ├── jump_compound_results_pattern_aug_070624.csv
│   │   └── taorf_results_pattern_aug_070624.csv
│   │
│   ├── virtual_screen_results_202407/
│   │   ├── lincs_screen_results.xlsx    [2.2 → Manuscript]
│   │   ├── CDRP_screen_results.xlsx
│   │   ├── jump_orf_screen_results.xlsx
│   │   ├── jump_crispr_screen_results.xlsx
│   │   ├── jump_compound_screen_results.xlsx
│   │   └── taorf_screen_results.xlsx
│   │
│   ├── Figure3.pdf                      [1.0 → Manuscript]
│   ├── SuppFigure2.pdf
│   ├── Figure4b.pdf
│   ├── Figure4c.pdf
│   └── dendrogram_mito.pdf
│
└── per_site_aggregated_profiles_newpattern_2/
    ├── lincs/                           [2.0 Step 3 → Step 4]
    │   ├── {batch1}_site_agg_profiles.csv.gz
    │   └── {batch2}_site_agg_profiles.csv.gz
    ├── CDRP/
    ├── jump_orf/
    ├── jump_crispr/
    ├── jump_compound/
    └── taorf/
```

---

## Critical Algorithm: MITO-SLOPE Calculation

This algorithm is central to both patient analysis (1.0) and virtual screening (2.0 Step 4).

```python
def find_end_slope2(data, height=None):
    """
    Calculate MITO-SLOPE: slope of mitochondrial radial distribution at cell edge

    Args:
        data: Array of radial distribution values (12 bins from nucleus edge to cell edge)
              Values are control-subtracted and z-scored

    Returns:
        tuple: (last_peak_ind, slope)
            last_peak_ind: Index of last extremum
            slope: Change in z-score per bin from last extremum to cell edge
    """

    # 1. SMOOTH DATA (reduce noise)
    data = savgol_filter(data, window_length=5, polyorder=3)

    # 2. FIND EXTREMA
    min_max_indc = [np.argmax(data), np.argmin(data)]

    # 3. FILTER TO VALID INDICES (exclude last 2 bins - edge artifacts)
    last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]

    IF last_peak_ind0 == []:
        RETURN (0, 0)  # No valid extremum

    # 4. SELECT LAST (RIGHTMOST) EXTREMUM
    last_peak_ind = max(last_peak_ind0)

    # 5. CALCULATE SLOPE FROM LAST EXTREMUM TO EDGE
    slope = (data[-1] - data[last_peak_ind]) / (len(data) - last_peak_ind - 1)

    RETURN (last_peak_ind, slope)

INTERPRETATION:
- slope > 0: Mitochondria enriched at cell periphery (more dispersion)
- slope < 0: Mitochondria depleted at cell periphery (perinuclear clustering)
- slope ≈ 0: Relatively uniform distribution

PATIENT RESULTS:
- Psychosis (BP+SZ+SZA): slope < 0 (p = 0.0069) ← Perinuclear clustering
- MDD: slope > 0 (p = 0.18, NS) ← Trend toward peripheral dispersion
- Controls: slope ≈ 0 (baseline)

SCREENING GOAL:
Find compounds/genes that:
- Increase slope (positive hits): Could treat psychosis
- Decrease slope (negative hits): Could treat depression
```

---

## Dataset-Specific Parameters Reference

```python
DATASET_CONFIGS = {
    'lincs': {
        'cell_line': 'A549',
        'n_compounds': 1571,
        'n_doses': 6,
        'replicates': 5,
        'pert_col': 'Metadata_pert_id_dose',
        'key_col': 'Metadata_pert_name',
        'controls': 'Metadata_pert_type == "control"',
        'bh_target_critical': 0.00761,
        'bh_orth_critical': 0.04829,
        'orth_features': 7
    },
    'CDRP': {
        'cell_line': 'U2OS',
        'n_compounds': 30430,
        'n_doses': 1,
        'replicates': 4,
        'pert_col': 'Metadata_Sample_Dose',
        'key_col': 'Metadata_pert_id',
        'controls': 'Metadata_pert_type == "control"',
        'bh_target_critical': 0.00924,
        'bh_orth_critical': 0.04812,
        'orth_features': 23
    },
    'jump_orf': {
        'cell_line': 'U2OS',
        'n_genes': 12609,
        'pert_col': 'Metadata_JCP2022',
        'key_col': 'Metadata_Symbol',
        'controls': 'Metadata_Symbol in ["LacZ", "BFP", "HcRed", "LUCIFERASE"]',
        'bh_target_critical': 0.00198,
        'bh_orth_critical': 0.03752,
        'orth_features': 40
    },
    'jump_crispr': {
        'cell_line': 'U2OS',
        'n_genes': 7975,
        'pert_col': 'Metadata_JCP2022',
        'key_col': 'Metadata_Symbol',
        'controls': 'Metadata_Symbol == "non-targeting"',
        'bh_target_critical': 0.007,
        'bh_orth_critical': 0.04785,
        'orth_features': 14
    },
    'jump_compound': {
        'cell_line': 'U2OS',
        'n_compounds': 116000,
        'pert_col': 'Metadata_JCP2022',
        'key_col': 'Metadata_InChIKey',
        'controls': 'DMSO',
        'bh_target_critical': 0.00365,
        'bh_orth_critical': 0.04746,
        'orth_features': 14
    },
    'taorf': {
        'cell_line': 'U2OS',
        'n_genes': 323,
        'replicates': 5,
        'pert_col': 'Metadata_broad_sample',
        'key_col': 'Metadata_gene_name',
        'controls': 'Metadata_gene_name in ["LacZ", "Luciferase"]',
        'bh_target_critical': 0.00409,
        'bh_orth_critical': 0.03287,
        'orth_features': 14
    }
}
```

---

## Common Pitfalls & Design Decisions

### Why Per-Plate Statistics?

**Problem:** Batch effects between plates can dominate biological signal

**Solution:** Calculate statistics independently per plate, then aggregate
- Standardize within plate (z-score)
- Test perturbation vs control within same plate
- Select median plate for final values

### Why Cohen's d as Primary Metric?

**Problem:** P-values depend on sample size, making comparisons across datasets difficult

**Solution:** Use Cohen's d (effect size) as primary ranking metric
- Independent of sample size
- Directly interpretable (standardized mean difference)
- Allows comparison across datasets

### Why Orthogonal Filter?

**Problem:** Compounds affecting mitochondria may be toxic and affect everything

**Solution:** Filter to perturbations that specifically affect MITO-SLOPE
- Target filter: p_slope_std < BH_critical (significant change in slope)
- Orthogonal filter: p_orth_std > BH_critical (no significant change in other features)
- Result: Perturbations with specific mitochondrial localization effect

### Why Two Standardization Steps in 2.0 Step 4?

1. **First standardization (line 1134):** Per-plate z-score of raw features
   - Makes features comparable within plate
   - Removes plate-level batch effects

2. **Control subtraction (line 1141):** Subtract control mean per plate
   - Removes plate-specific baseline
   - Creates control-relative pattern

3. **Slope calculation (line 1156):** Apply find_end_slope2
   - Converts 12-bin pattern to single slope value

4. **Second standardization (line 1167):** Per-plate z-score including slope
   - Makes slope comparable across plates
   - Enables cross-plate aggregation

---

## Reproducibility Checklist

To reproduce this analysis from scratch:

### Prerequisites
```bash
# Python environment
conda create -n mito_analysis python=3.9
conda activate mito_analysis
pip install -r requirements.txt

# Key packages
- pandas
- numpy
- scipy
- scikit-learn
- matplotlib
- seaborn
- blitzgsea
- singlecell-morph (custom package)
```

### Data Access
```
1. Patient fibroblast data:
   - Contact: mhaghigh@broadinstitute.org
   - Location: S3 bucket (see download_data.sh)
   - Size: 4.43 GB

2. Public Cell Painting data:
   - LINCS: Cell Painting Gallery cpg0004
   - CDRP: Cell Painting Gallery cpg0012
   - JUMP: Cell Painting Gallery cpg0016
   - TA-ORF: Cell Painting Gallery cpg0017
```

### Execution Order
```
1. Run download_data.sh (if data not locally available)
2. Run notebooks/1.0-mh-feat-importance.py
3. Run notebooks/2.0-mh-virtual-screen.py (Steps 1-4 sequentially)
4. Run notebooks/2.1-mh-set-enrichment-analysis.py
5. Run notebooks/2.2-mh-check-vs-lists.py
```

### Expected Runtime
```
- 1.0: ~30 minutes (depends on patient data size)
- 2.0 Step 1: ~5 minutes (metadata only)
- 2.0 Step 2: ~2 hours per dataset (depends on backend size)
- 2.0 Step 3: ~4 hours per dataset (SQLite I/O intensive)
- 2.0 Step 4: ~30 minutes per dataset
- 2.1: ~15 minutes (API calls to Enrichr)
- 2.2: ~5 minutes
```

---

## Version Control & Data Provenance

```
CODE VERSION: Git commit from 2025_Haghighi_Mito repository

DATA VERSIONS:
- Patient fibroblasts: patient_labels_updatedSept302025.xlsx
- LINCS: cpg0004 (2022 release)
- CDRP: cpg0012 (2017 release)
- JUMP: cpg0016 (2024 release)
- TA-ORF: cpg0017 (2017 release)

ALGORITHM VERSIONS:
- MITO-SLOPE: Version defined in 1.0 (lines 788-840)
- CellProfiler features: Extracted via CellProfiler 3/4
- Statistical tests: scipy.stats (2-sample t-test, Hotelling's T²)
- Enrichment: blitzGSEA (version from requirements)
- BH correction: Custom implementation (lines 266-278 in 2.2)
```

---

## Contact & Support

For questions about:
- **Patient data:** Marzieh Haghighi (mhaghigh@broadinstitute.org)
- **Analysis pipeline:** Shantanu Singh (shsingh@broadinstitute.org)
- **Cell Painting data:** Beth Cimini (bcimini@broadinstitute.org)

---

*Document generated: 2025-10-17*
*Pipeline version: 2024_Haghighi_Mito commit 786587a*
