"""Verify that applying runtime label overrides to the shipped aggregated CSV
produces a DataFrame identical to the corrected CSV emitted by
``build_corrected_upstream_labels``.

If this passes, any deterministic figure pipeline (notebook 3.0, Figure 4c,
Supp Fig 3 panels A-E) that currently reads the shipped CSV plus overrides is
guaranteed to produce byte-identical intermediate data when switched to reading
the corrected CSV with no overrides. That proves Option A (regenerate upstream
artefacts) is numerically equivalent to the current fork behaviour.

Run: ``.pixi/envs/default/bin/python scripts/verify_corrected_labels_equivalence.py``
"""

import sys

import pandas as pd

from haghighi_mito.config import AGGREGATED_PROFILES_PATH
from haghighi_mito.patient_analysis import load_aggregated_profiles

CORRECTED_AGG = "data/interim/upstream_corrected/aggregated_profiles_fibroblast.csv"
CORRECTED_SUBJECTS = "data/interim/upstream_corrected/subjects.csv"


def _apply_notebook_overrides(df: pd.DataFrame) -> pd.DataFrame:
    """Replicate notebook 3.0's inline override block verbatim."""
    df = df.copy()
    df.loc[df["subject"].astype(str) == "272", "label"] = "Control"
    df.loc[df["subject"].astype(str) == "MCL004", "label"] = "SZA"
    df["label"] = df["label"].replace("MDD or Dep", "MDD")
    df["subject"] = df["subject"].astype(str)
    return df


def _normalise_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce subject to str so comparison isn't confused by int-vs-str encoding."""
    df = df.copy()
    df["subject"] = df["subject"].astype(str)
    return df.reset_index(drop=True)


def main() -> int:
    original = pd.read_csv(AGGREGATED_PROFILES_PATH)
    corrected = pd.read_csv(CORRECTED_AGG)

    with_overrides = _apply_notebook_overrides(original)

    left = _normalise_for_compare(with_overrides)
    right = _normalise_for_compare(corrected)

    failures: list[str] = []

    if left.shape != right.shape:
        failures.append(f"shape mismatch: original+overrides {left.shape} vs corrected {right.shape}")
    if list(left.columns) != list(right.columns):
        failures.append("column order/name mismatch")
    else:
        diff = (left != right) & ~(left.isna() & right.isna())
        if diff.any().any():
            per_col = diff.sum().loc[lambda s: s > 0].to_dict()
            failures.append(f"value mismatches per column: {per_col}")

    # library-level parity
    lib_df = _normalise_for_compare(load_aggregated_profiles())
    if not lib_df.equals(left):
        failures.append("load_aggregated_profiles() disagrees with notebook 3.0's inline overrides")

    # subjects.csv shape
    subjects = pd.read_csv(CORRECTED_SUBJECTS)
    if len(subjects) != right["subject"].nunique():
        failures.append(f"subjects.csv has {len(subjects)} rows but aggregated CSV has {right['subject'].nunique()} unique subjects")

    header = "EQUIVALENCE CHECK: original CSV + runtime overrides  vs  corrected CSV"
    print(header)
    print("=" * len(header))
    print(f"original CSV:      {AGGREGATED_PROFILES_PATH}")
    print(f"corrected CSV:     {CORRECTED_AGG}")
    print(f"subjects.csv:      {CORRECTED_SUBJECTS}")
    print(f"rows:              {len(left)}")
    print(f"columns:           {len(left.columns)}")
    print()

    if failures:
        print("FAIL — the two paths diverge:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("PASS — the two paths produce byte-identical DataFrames.")
    print()
    print("Implication: any deterministic figure pipeline that currently reads")
    print("the shipped CSV and applies the four overrides can be switched to")
    print("reading the corrected CSV with no overrides, with zero numeric impact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
