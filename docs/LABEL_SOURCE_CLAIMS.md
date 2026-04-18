# Patient-label source claims — verification sheet

**Purpose.** This document enumerates every factual claim made about how patient
disease labels (`Control`, `BP`, `SZ`, `SZA`, `MDD`) flow through three
pipelines: **Erin's** supplemental cellsize notebook, **Marzieh's**
phenotype-discovery notebooks, and **our fork**. Each claim is paired with an
exact file path, an expected string/pattern, and a one-line shell/python check
that a second agent can run to confirm or refute it.

**Scope.** Labels only. Feature computation, QC, slope, statistics, and figures
are out of scope.

**How to verify.** Run each `Check` block from the repo root
(`/Users/shsingh/Documents/GitHub/mito/2025_Haghighi_Mito`). A claim is
**verified** iff its check produces the `Expected` output exactly. If a check
fails, the claim is wrong — flag it, do not patch the check.

**Verifier posture.** Do not trust this document. Re-run every check. If a file
path has moved or a line number has drifted, report the drift rather than
silently adjusting.

---

## Environment assumptions

- Repo root: `/Users/shsingh/Documents/GitHub/mito/2025_Haghighi_Mito`
- Upstream worktree exists at `.worktrees/upstream/` with `remotes/upstream/HEAD → main`
- Erin's notebook lives on upstream default branch at commit `ef8f26b`
  (`update singlecell plots`). It is read via `git show` because the current
  checkout of the worktree is an older branch (`upstream-main`) that predates
  the commit.
- Python on PATH must be able to `import json` (stdlib only for all checks).

---

## Section 1 — Erin's notebook (supplemental cellsize)

**Path on upstream main:**
`results/supplemental_single_cell_variability_cellsize/plot_singlecell_cellsize.ipynb`
at commit `ef8f26b`.

### Claim E1 — Erin's only label source is `aggregated_profiles_fibroblast.csv`

Erin reads labels from exactly one file: `aggregated_profiles_fibroblast.csv`.
She does **not** read any `.pkl`, she does **not** read any
`patient_labels*.xlsx`, she does **not** read `subjects.csv`.

**Check:**
```bash
git -C .worktrees/upstream show ef8f26b:results/supplemental_single_cell_variability_cellsize/plot_singlecell_cellsize.ipynb \
  | python3 -c "
import json, sys
nb = json.load(sys.stdin)
reads = []
for c in nb['cells']:
    if c['cell_type'] != 'code': continue
    src = ''.join(c['source'])
    for needle in ['.pkl', '.xlsx', 'subjects.csv', 'patient_labels', 'aggregated_profiles_fibroblast']:
        if needle in src:
            reads.append(needle)
print(sorted(set(reads)))
"
```
**Expected:** `['aggregated_profiles_fibroblast']`

### Claim E2 — Erin applies no subject-specific disease-label overrides

Erin does not apply the subject-specific disease-label corrections that appear
in Marzieh's notebooks and our fork: no `272→Control`, no `MCL004→SZA`, no
`370{E,F,H}→370`, and no `"MDD or Dep"→MDD` rewrite. She *does* manipulate the
`label` column later for plotting/category order and to append the composite
`psychosis` group; those are not patient-label source overrides.

**Check:**
```bash
git -C .worktrees/upstream show ef8f26b:results/supplemental_single_cell_variability_cellsize/plot_singlecell_cellsize.ipynb \
  | python3 -c "
import json, sys
nb = json.load(sys.stdin)
hits = []
for c in nb['cells']:
    if c['cell_type'] != 'code': continue
    src = ''.join(c['source'])
    for needle in [
        \"replace(['370E','370F','370H'],'370')\",
        \"replace('MDD or Dep','MDD')\",
        \"loc[df_1['subject']=='272','label']='Control'\",
        \"loc[df_1['subject']=='MCL004','label']='SZA'\",
    ]:
        if needle in src:
            hits.append(needle)
print(sorted(set(hits)))
"
```
**Expected:** `[]`

### Claim E3 — Erin's label join is `df_new.merge(ag_df[['subject','label']], on='subject', how='left')`

The merge is a left join of `(subject,label)` from the aggregated CSV onto the
SQL-loaded cell table. No other join carries labels.

**Check:**
```bash
git -C .worktrees/upstream show ef8f26b:results/supplemental_single_cell_variability_cellsize/plot_singlecell_cellsize.ipynb \
  | python3 -c "
import json, sys
nb = json.load(sys.stdin)
for c in nb['cells']:
    if c['cell_type'] != 'code': continue
    src = ''.join(c['source'])
    if \"ag_df[['subject','label']]\" in src:
        print('FOUND')
        break
else:
    print('MISSING')
"
```
**Expected:** `FOUND`

### Claim E4 — Erin excludes subject `'370'` from the SQL loop

Before looping SQL backends, Erin builds `population` and drops `'370'` with
the comment "drop 370 because it was split into 370A-370H".

**Check:**
```bash
git -C .worktrees/upstream show ef8f26b:results/supplemental_single_cell_variability_cellsize/plot_singlecell_cellsize.ipynb \
  | grep -c "x for x in population if x != '370'"
```
**Expected:** `1`

---

## Section 2 — Marzieh's notebooks (phenotype discovery)

**Paths (current upstream worktree HEAD, branch `upstream-main`):**
- `.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb`
- `.worktrees/upstream/phenotype_discovery/1.feature_inspection.ipynb`
- `.worktrees/upstream/phenotype_discovery/subjects.csv`

### Claim M1 — Marzieh's cell 4 reads four sources, in this order

In `2.slope_analysis.ipynb` cell 4, these reads occur (in this order):

1. `pd.read_pickle(dataDir+'single_cell_with_annot.pkl', ...)` → `df_1_0_0`
2. Per-subject loop `readSingleCellData_sqlalch(.../<si>.sqlite)` → concat → `df_new`
3. `pd.read_pickle(dataDir+'single_cell_with_annot_allFeatures.pkl', ...)` → `df_1_0`
4. `pd.read_excel(...patient_labels_updatedSept302025.xlsx)` → `disease_labels`

**Check:**
```bash
python3 -c "
import json
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src = ''.join(nb['cells'][4]['source'])
ordered_needles = [
    \"read_pickle(dataDir+'single_cell_with_annot.pkl'\",
    'readSingleCellData_sqlalch',
    \"read_pickle(dataDir+'single_cell_with_annot_allFeatures.pkl'\",
    'patient_labels_updatedSept302025.xlsx',
]
positions = [src.find(n) for n in ordered_needles]
print('all_found:', all(p >= 0 for p in positions))
print('in_order :', positions == sorted(positions))
"
```
**Expected:**
```
all_found: True
in_order : True
```

### Claim M2 — Marzieh's label merges are pickle-first, then xlsx

The label column is first attached via
`pd.merge(df_1, df_1_0[['subject','label']].drop_duplicates(), on='subject', how='left')`
and then `disease_labels` is left-joined on `subject`.

**Check:**
```bash
python3 -c "
import json
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src = ''.join(nb['cells'][4]['source'])
m1 = src.find(\"df_1_0[['subject','label']].drop_duplicates()\")
m2 = src.find('pd.merge(df_1,disease_labels,on=[\"subject\"]')
print('both_present:', m1 > 0 and m2 > 0)
print('pickle_then_xlsx:', m1 < m2)
"
```
**Expected:**
```
both_present: True
pickle_then_xlsx: True
```

### Claim M3 — Marzieh applies four inline overrides in cell 4

After the two label merges, cell 4 executes:
- `df_1['subject'] = df_1['subject'].replace(['370E','370F','370H'],'370')`
- `df_1.loc[df_1['subject']=='272','label']='Control'`
- `df_1.loc[df_1['subject']=='MCL004','label']='SZA'`
- `df_1['label'] = df_1['label'].replace('MDD or Dep','MDD')` (later in the same cell, before scaling)

**Check:**
```bash
python3 -c "
import json
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src = ''.join(nb['cells'][4]['source'])
for needle in [
    \"replace(['370E','370F','370H'],'370')\",
    \"df_1.loc[df_1['subject']=='272','label']='Control'\",
    \"df_1.loc[df_1['subject']=='MCL004','label']='SZA'\",
    \"replace('MDD or Dep','MDD')\",
]:
    print(f'{needle!r:70s}', '->', needle in src)
"
```
**Expected:** four lines, each ending with `-> True`.

### Claim M4 — The `to_csv(...aggregated_profiles_fibroblast.csv)` call in cell 4 is guarded by `if 0:`

Marzieh has the write, but by default it does **not** execute on notebook
re-run. The CSV on S3 was produced at some earlier moment when the guard was
removed (or run manually).

This is a correction to an earlier claim in-session that said Marzieh writes
the CSV "at the end of cell 4." She *has* the write but it is gated off.

**Check:**
```bash
python3 -c "
import json, re
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src = ''.join(nb['cells'][4]['source'])
# Find any line matching the aggregated-profiles to_csv
for line in src.splitlines():
    if 'aggregated_profiles_fibroblast.csv' in line and 'to_csv' in line:
        print(repr(line))
# and show a small surrounding block
idx = src.find('aggregated_profiles_fibroblast.csv')
print('---context---')
print(src[max(0,idx-80):idx+120])
"
```
**Expected:** the printed line is `"    df_1_avg_persub.to_csv(dataDir+'/aggregated_profiles_fibroblast.csv',index=False)"` and the context block shows the preceding `if 0:` guard.

### Claim M5 — `subjects.csv` is generated by a commented-out line in cell 6

In `2.slope_analysis.ipynb` cell 6 there is a single line
`# df_1_avg_persub[['subject','label']].to_csv('subjects.csv',index=False)` —
commented out. `phenotype_discovery/subjects.csv` therefore has to have been
produced by uncommenting and running it once, not by any automated step.

**Check:**
```bash
python3 -c "
import json
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src6 = ''.join(nb['cells'][6]['source'])
for line in src6.splitlines():
    if 'subjects.csv' in line:
        print(repr(line))
"
```
**Expected:** exactly one line, beginning with `'# ` (a comment).

### Claim M6 — Marzieh's cell 26 reads an older labels file

Cell 26 of `2.slope_analysis.ipynb` contains
`pd.read_excel(.../patient_labels.xlsx)` (no `_updatedSept302025`). This is a
separate, later merge onto `data_phs0`, used only to pull `D1` (a finer
diagnostic subtype column) for the boxplot in that cell. It does **not**
overwrite the `label` column. The earlier check for this claim was too loose:
it treated a commented read of `data_phs['label']` as if it were an
assignment.

**Check:**
```bash
python3 -c "
import json, re
nb = json.load(open('.worktrees/upstream/phenotype_discovery/2.slope_analysis.ipynb'))
src = ''.join(nb['cells'][26]['source'])
print('reads_old_xlsx   :', 'patient_labels.xlsx' in src and 'updatedSept302025' not in src.split('patient_labels.xlsx')[0].rsplit('pd.read_excel',1)[-1])
print('merges_on_subject:', 'on=[\"subject\"]' in src)
print('writes_label     :', bool(re.search(r\"\\[['\\\"]label['\\\"]\\]\\s*=(?!=)\", src)))
"
```
**Expected:**
```
reads_old_xlsx   : True
merges_on_subject: True
writes_label     : False
```
(The third check now asks specifically whether cell 26 assigns to the
`'label'` column — it should not. If it does, claim M6 is wrong and needs
updating.)

### Claim M7 — `1.feature_inspection.ipynb` repeats the same cell-4 pattern

Marzieh's feature-inspection notebook has a cell with the same four reads
and the same four overrides. The cell index may differ from `2.slope_analysis`;
the check below looks across all code cells.

**Check:**
```bash
python3 -c "
import json
nb = json.load(open('.worktrees/upstream/phenotype_discovery/1.feature_inspection.ipynb'))
matches = 0
for c in nb['cells']:
    if c['cell_type'] != 'code': continue
    src = ''.join(c['source'])
    if (\"single_cell_with_annot.pkl\" in src
        and \"single_cell_with_annot_allFeatures.pkl\" in src
        and \"patient_labels_updatedSept302025.xlsx\" in src
        and \"df_1.loc[df_1['subject']=='272','label']='Control'\" in src
        and \"df_1.loc[df_1['subject']=='MCL004','label']='SZA'\" in src):
        matches += 1
print('cells_with_full_pattern:', matches)
"
```
**Expected:** `cells_with_full_pattern: 1`

---

## Section 3 — Our fork

**Paths (all relative to repo root):**
- `haghighi_mito/config.py`
- `haghighi_mito/patient_analysis.py`
- `notebooks/3.0-mh-slope-analysis.py`
- `notebooks/1.0-mh-feat-importance.py`
- Downloaded data: `data/external/tables/curated_2025-10-25/` and
  `data/external/mito_project/workspace/singleCellData/`

### Claim F1 — Two label CSVs are tracked/downloadable

`haghighi_mito/config.py` defines:
- `PATIENT_LABELS_PATH = EXTERNAL_TABLES_DIR / "patient_labels_updatedSept302025.csv"`
- `PATIENT_LABELS_FULL_PATH = EXTERNAL_TABLES_DIR / "patient_labels_updatedSept302025_full.csv"`

Both files currently exist on disk under
`data/external/tables/curated_2025-10-25/`.

**Check:**
```bash
grep -nE "PATIENT_LABELS_(FULL_)?PATH" haghighi_mito/config.py
ls data/external/tables/curated_2025-10-25/patient_labels_updatedSept302025*.csv
```
**Expected:** `PATIENT_LABELS_PATH` and `PATIENT_LABELS_FULL_PATH` assignments
to the two filenames, and `ls` shows both `.csv` files.

### Claim F2 — `load_patient_labels()` prefers the `_full.csv` when present

`patient_analysis.py::load_patient_labels` uses
`PATIENT_LABELS_FULL_PATH if PATIENT_LABELS_FULL_PATH.exists() else PATIENT_LABELS_PATH`,
renames `ID → subject`, and stringifies `subject`.

**Check:**
```bash
python3 -c "
import re
src = open('haghighi_mito/patient_analysis.py').read()
m = re.search(r'def load_patient_labels.*?(?=\ndef |\Z)', src, re.S)
assert m, 'function not found'
body = m.group(0)
print('prefers_full :', 'PATIENT_LABELS_FULL_PATH.exists()' in body)
print('renames_id   :', 'rename(columns={\"ID\": \"subject\"})' in body)
print('stringifies  :', 'astype(str)' in body)
"
```
**Expected:**
```
prefers_full : True
renames_id   : True
stringifies  : True
```

### Claim F3 — `load_single_cell_data()` replicates Marzieh's label cascade

In `patient_analysis.py::load_single_cell_data`:
1. Reads all SQLite backends and concatenates → `df_new`
2. Reads `single_cell_with_annot_allFeatures.pkl` → `df_allfeatures`
3. Merges label via `df_allfeatures[['subject','label']].drop_duplicates()`
4. Merges disease labels via `load_patient_labels()`
5. Applies `replace(['370E','370F','370H'],'370')`
6. Applies `df_raw.loc[df_raw['subject']=='272','label']='Control'`
7. Applies `df_raw.loc[df_raw['subject']=='MCL004','label']='SZA'`

Note: step 7 is present but **`"MDD or Dep" → "MDD"` is NOT applied in
`load_single_cell_data`** — it is applied only later in `preprocess_single_cells`
(line with `df["label"] = df["label"].replace("MDD or Dep", "MDD")`). This
asymmetry is worth flagging.

**Check:**
```bash
python3 -c "
import re
src = open('haghighi_mito/patient_analysis.py').read()
m = re.search(r'def load_single_cell_data.*?(?=\ndef |\Z)', src, re.S)
body = m.group(0)
ordered = [
    'readSingleCellData_sqlalch',
    'single_cell_with_annot_allFeatures.pkl',
    \"df_allfeatures[['subject', 'label']].drop_duplicates()\",
    'load_patient_labels()',
    \"replace(['370E', '370F', '370H'], '370')\",
    'df_raw[\"subject\"] == \"272\"',
    'df_raw[\"subject\"] == \"MCL004\"',
]
positions = [body.find(n) for n in ordered]
for n, p in zip(ordered, positions):
    print(f'{p:5d}  {n}')
print('all_found:', all(p >= 0 for p in positions))
print('in_order :', positions == sorted(positions))
# And confirm MDD-or-Dep is NOT in load_single_cell_data (it lives in preprocess)
print('mdd_replace_in_load_single :', 'MDD or Dep' in body)
pp = re.search(r'def preprocess_single_cells.*?(?=\ndef |\Z)', src, re.S).group(0)
print('mdd_replace_in_preprocess  :', 'MDD or Dep' in pp)
"
```
**Expected:**
- `all_found: True`
- `in_order : True`
- `mdd_replace_in_load_single : False`
- `mdd_replace_in_preprocess  : True`

### Claim F4 — `load_aggregated_profiles()` reads the aggregated CSV and applies three idempotent overrides

`patient_analysis.py::load_aggregated_profiles` reads
`AGGREGATED_PROFILES_PATH`, applies `272→Control`, `MCL004→SZA`,
`"MDD or Dep"→MDD`, and does **not** collapse `370E/F/H` (they are not present
in the aggregated CSV — already collapsed upstream).

**Check:**
```bash
python3 -c "
import re
src = open('haghighi_mito/patient_analysis.py').read()
m = re.search(r'def load_aggregated_profiles.*?(?=\ndef |\Z)', src, re.S)
body = m.group(0)
print('reads_csv     :', 'pd.read_csv(AGGREGATED_PROFILES_PATH)' in body)
print('override_272  :', '\"272\"' in body and \"'Control'\" in body)
print('override_MCL  :', '\"MCL004\"' in body and \"'SZA'\" in body)
print('mdd_replace   :', \"'MDD or Dep'\" in body and \"'MDD'\" in body)
print('no_370_collapse:', '370' not in body)
"
```
**Expected:**
```
reads_csv     : True
override_272  : True
override_MCL  : True
mdd_replace   : True
no_370_collapse: True
```

### Claim F5 — Notebook 3.0 reads labels from aggregated CSV (not pickle, not xlsx)

`notebooks/3.0-mh-slope-analysis.py` loads labels only from
`AGGREGATED_PROFILES_PATH`, applies the same three overrides as
`load_aggregated_profiles`, and uses those labels both for the aggregated-profile
analysis and as the label join onto the SQL single-cell parquet (panels A–E).

**Check:**
```bash
grep -nE "read_csv|read_pickle|read_excel|AGGREGATED_PROFILES_PATH|FIBROBLAST_SINGLECELL_PARQUET" \
  notebooks/3.0-mh-slope-analysis.py
```
**Expected:** `read_csv(AGGREGATED_PROFILES_PATH)` appears; no `read_pickle`,
no `read_excel` (except possibly inside comments). A `read_parquet` on
`FIBROBLAST_SINGLECELL_PARQUET` is expected for the raw-micron panels but
carries **no** label column — labels are merged in from the aggregated CSV.

### Claim F6 — Notebook 1.0 uses the heavy `load_single_cell_data()` path

`notebooks/1.0-mh-feat-importance.py` imports `load_single_cell_data` and
calls it. It does not call `load_aggregated_profiles`.

**Check:**
```bash
grep -nE "load_single_cell_data|load_aggregated_profiles" notebooks/1.0-mh-feat-importance.py
```
**Expected:** one import line and one `load_single_cell_data()` call; no
reference to `load_aggregated_profiles`.

---

## Section 4 — The cross-pipeline claim

### Claim X1 — The shipped aggregated CSV only partially reflects Marzieh's four label corrections

Marzieh's cell 4 applies all four corrections explicitly. Our fork also applies
all four corrections, split across `load_single_cell_data`,
`preprocess_single_cells`, `load_aggregated_profiles`, and notebook 3.0.
However, the shipped `aggregated_profiles_fibroblast.csv` on disk only
reflects two of those four corrections: it has collapsed `370` and no
`"MDD or Dep"` label, but it still carries `272→SZA` and `MCL004→MDD`.
Therefore Erin's notebook, which reads labels only from that CSV and applies no
subject-specific rewrites, does **not** end up with the same final labels as
Marzieh / our fork for those two subjects.

| Correction | Marzieh (cell 4) | Erin | Our fork (heavy) | Our fork (light) |
|---|---|---|---|---|
| `370{E,F,H} → 370` | explicit | inherited via aggregated CSV | explicit (`load_single_cell_data`) | inherited via aggregated CSV |
| `272 → Control`    | explicit | **not present in shipped CSV** | explicit | explicit (idempotent) |
| `MCL004 → SZA`     | explicit | **not present in shipped CSV** | explicit | explicit (idempotent) |
| `"MDD or Dep" → MDD` | explicit | inherited via aggregated CSV | explicit in `preprocess_single_cells` | explicit (idempotent) |

This supersedes the earlier hypothesis that the aggregated CSV on disk had
already inherited all four of Marzieh's corrections. The local file refutes
that hypothesis.

**Check (data-dependent — requires the downloaded CSV):**
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('data/external/mito_project/workspace/singleCellData/aggregated_profiles_fibroblast.csv')
print('rows                :', len(df))
print('unique subjects     :', df['subject'].nunique())
print('label values        :', sorted(df['label'].dropna().unique().tolist()))
print('has_370E_F_H        :', any(s in set(df['subject'].astype(str).unique()) for s in ['370E','370F','370H']))
print('has_collapsed_370   :', '370' in set(df['subject'].astype(str).unique()))
print('272_label           :', df.loc[df['subject'].astype(str)=='272','label'].unique().tolist())
print('MCL004_label        :', df.loc[df['subject'].astype(str)=='MCL004','label'].unique().tolist())
print('has_MDD_or_Dep      :', 'MDD or Dep' in set(df['label'].dropna().astype(str).unique()))
"
```
**Expected:**
- `has_370E_F_H        : False`
- `has_collapsed_370   : True`
- `272_label           : ['SZA']`
- `MCL004_label        : ['MDD']`
- `has_MDD_or_Dep      : False`
- `label values`: subset of `['BP','Control','MDD','SZ','SZA']`

If any of these fail, claim X1 is wrong for the CSV on disk and needs updating.

---

## Section 5 — Open questions (explicitly NOT yet claims)

These are things **not yet verified** — they are the targets for the next
session's work. They are listed here so a verifier can flag any claim that
drifts into these areas without evidence.

1. **Goal 1: numerical match with Erin.** It is *plausible* that our light
   path (`load_aggregated_profiles` + notebook 3.0) matches Erin's
   per-subject means exactly, because both consume the same aggregated CSV
   and both join on `subject`. But this has not been checked by comparing
   per-subject slopes, per-subject raw-micron means, or any figure.
2. **Goal 2: numerical match with Marzieh.** It is *plausible* that our heavy
   path (`load_single_cell_data` + `preprocess_single_cells`) matches
   Marzieh's per-subject slopes exactly. This requires running notebook 1.0
   end-to-end and comparing to the `slope` column in
   `aggregated_profiles_fibroblast.csv`.
3. **`subjects.csv` coverage.** `phenotype_discovery/subjects.csv` has 168
   subjects while the pickle and SQL have 176. The 8-subject delta has not
   been enumerated or reconciled in this document.
4. **`single_cell_with_annot.pkl` vs `_allFeatures.pkl`.** Marzieh reads
   both. Only `_allFeatures.pkl` is used for the label merge. The role of
   the non-`allFeatures` pickle in cell 4 (assigned to `df_1_0_0` and never
   re-referenced in the cells inspected) is not established here.

A verifier should refuse to treat any of the above as true without an
explicit, checkable claim added to Sections 1–4.
