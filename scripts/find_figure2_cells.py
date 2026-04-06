#!/usr/bin/env python3
"""Identify the representative single cell within each Figure 2 image.

Usage: python scripts/find_figure2_cells.py <csv_path> <output_csv>

For each row in figure2_representative_cells.csv, loads only the Cells compartment
from the subject's SQLite, filters to the chosen image, computes per-cell slope on
raw MeanFrac values, and picks the cell closest to the within-image median slope.
Outputs ObjectNumber, center coordinates, and bounding box for cropping.

No QC filtering or standardisation is needed: all cells are in the same image so
raw values are directly comparable for ranking purposes.
"""
import sys

import sqlite3

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

from haghighi_mito.config import SQLITE_DATA_DIR

TARGET_COLUMNS = [
    f"Cells_RadialDistribution_MeanFrac_mito_tubeness_CEN_{i}of12"
    for i in range(1, 13)
]

LOCATION_COLS = [
    "ObjectNumber",
    "Cells_Location_Center_X",
    "Cells_Location_Center_Y",
    "Cells_AreaShape_BoundingBoxMinimum_X",
    "Cells_AreaShape_BoundingBoxMinimum_Y",
    "Cells_AreaShape_BoundingBoxMaximum_X",
    "Cells_AreaShape_BoundingBoxMaximum_Y",
]


def cell_slope(row):
    data = savgol_filter(row.values, window_length=5, polyorder=3)
    candidates = [i for i in [np.argmax(data), np.argmin(data)] if i < len(data) - 2]
    if not candidates:
        return np.nan
    peak = np.max(candidates)
    return (data[-1] - data[peak]) / (len(data) - peak - 1)


csv_path, output_csv = sys.argv[1], sys.argv[2]
rep_images = pd.read_csv(csv_path)

results = []
for _, row in rep_images.iterrows():
    subject, label, image_number = row["subject"], row["label"], int(row["ImageNumber"])
    print(f"Processing {subject} ({label}), image {image_number}...")

    sqlite_path = str(SQLITE_DATA_DIR / subject / f"{subject}.sqlite")
    with sqlite3.connect(sqlite_path) as con:
        df_img = pd.read_sql(
            "SELECT * FROM Cells WHERE ImageNumber = ?", con, params=(image_number,)
        )

    missing = [c for c in TARGET_COLUMNS if c not in df_img.columns]
    if df_img.empty or missing:
        print(f"  WARNING: image {image_number} not found or missing {len(missing)} columns — skipping")
        continue

    df_img["cell_slope"] = df_img[TARGET_COLUMNS].apply(cell_slope, axis=1)
    df_img = df_img.dropna(subset=["cell_slope"])

    img_median = df_img["cell_slope"].median()
    best = df_img.loc[(df_img["cell_slope"] - img_median).abs().idxmin()]

    result = {
        "subject": subject,
        "label": label,
        "ImageNumber": image_number,
        "FileName_Mito": row["FileName_Mito"],
        "image_median_slope": img_median,
        "cell_slope": best["cell_slope"],
        "n_cells_in_image": len(df_img),
    }
    for col in LOCATION_COLS:
        result[col] = best[col] if col in best.index else np.nan

    results.append(result)
    print(f"  → ObjectNumber {int(best['ObjectNumber'])}, "
          f"center ({best['Cells_Location_Center_X']:.0f}, {best['Cells_Location_Center_Y']:.0f}), "
          f"slope {best['cell_slope']:.4f} (image median {img_median:.4f})")

out = pd.DataFrame(results)
out.to_csv(output_csv, index=False)
print(f"\nSaved to {output_csv}")
print(out[["subject", "label", "ImageNumber", "ObjectNumber",
           "Cells_Location_Center_X", "Cells_Location_Center_Y", "cell_slope"]].to_string(index=False))
