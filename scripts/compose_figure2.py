#!/usr/bin/env -S uv run --with tifffile --with matplotlib --with numpy --with pandas python3
"""Compose Figure 2: representative single cells per diagnosis group.

Usage: python scripts/compose_figure2.py <cells_csv> <images_dir> <output_pdf>

Layout: 5 rows (Control, BP, SZ, SZA, MDD) × 4 columns
        (DNA, Mito, Actin, composite).
        Each panel shows a white bounding box around the selected cell
        and a white × at the cell center.
Canonical subjects matching upstream: MCL162, MCL128, MCL015, 263, 287.

Channels: c1=DNA (blue), c2=Actin (green), c3=Mito (red).
Normalization: 0.5/99.5 percentile for all channels.
"""
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile

# Upstream-canonical one subject per group (from PROGRESS.md 2026-04-05 entry)
CANONICAL = {
    "BP": "MCL128",
    "Control": "MCL162",
    "MDD": "287",
    "SZ": "MCL015",
    "SZA": "263",
}
GROUP_ORDER = ["Control", "BP", "SZ", "SZA", "MDD"]

CHANNEL_COLOR = {
    "c1": np.array([0.2, 0.4, 1.0]),   # DNA → blue
    "c2": np.array([0.2, 0.9, 0.2]),   # Actin → green
    "c3": np.array([1.0, 0.2, 0.2]),   # Mito → red
}
COL_ORDER = ["c1", "c3", "c2", "merge"]
COL_TITLES = ["DNA", "Mito", "Actin", "Composite"]

BBOX_PAD_PX = 30  # padding around cell bounding box for the white rectangle


def normalize(img, plo=0.5, phi=99.5):
    lo, hi = np.percentile(img, plo), np.percentile(img, phi)
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1)


def make_merge(channels):
    h, w = next(iter(channels.values())).shape
    rgb = np.zeros((h, w, 3), dtype=float)
    for ch, img_norm in channels.items():
        rgb += img_norm[:, :, None] * CHANNEL_COLOR[ch]
    return np.clip(rgb, 0, 1)


def load_channels(images_dir, filename_mito):
    """Load all three channels, returning normalized arrays keyed by channel."""
    result = {}
    for ch in ("c1", "c2", "c3"):
        fname = filename_mito.replace("_c3_", f"_{ch}_")
        raw = tifffile.imread(str(images_dir / fname))
        result[ch] = normalize(raw)
    return result


cells_csv, images_dir, output_pdf = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
df = pd.read_csv(cells_csv)
df["subject"] = df["subject"].astype(str)

fig, axes = plt.subplots(
    len(GROUP_ORDER), len(COL_ORDER),
    figsize=(3 * len(COL_ORDER), 3 * len(GROUP_ORDER)),
    squeeze=False,
)

for row_i, group in enumerate(GROUP_ORDER):
    subject = CANONICAL[group]
    row_data = df[(df["subject"] == subject) & (df["label"] == group)]
    if row_data.empty:
        print(f"WARNING: no data for {subject} ({group}) — skipping row", file=sys.stderr)
        for ax in axes[row_i]:
            ax.axis("off")
        continue

    r = row_data.iloc[0]
    channels = load_channels(images_dir, r["FileName_Mito"])
    merge = make_merge(channels)

    cx = float(r["Cells_Location_Center_X"])
    cy = float(r["Cells_Location_Center_Y"])
    bx_min = int(r["Cells_AreaShape_BoundingBoxMinimum_X"])
    by_min = int(r["Cells_AreaShape_BoundingBoxMinimum_Y"])
    bx_max = int(r["Cells_AreaShape_BoundingBoxMaximum_X"])
    by_max = int(r["Cells_AreaShape_BoundingBoxMaximum_Y"])

    for col_i, col_key in enumerate(COL_ORDER):
        ax = axes[row_i][col_i]
        ax.axis("off")

        if col_key == "merge":
            ax.imshow(merge, interpolation="nearest", aspect="auto")
        else:
            ax.imshow(channels[col_key], cmap="gray", interpolation="nearest", aspect="auto")

        # White bounding box with padding
        rect_x = bx_min - BBOX_PAD_PX
        rect_y = by_min - BBOX_PAD_PX
        rect_w = (bx_max - bx_min) + 2 * BBOX_PAD_PX
        rect_h = (by_max - by_min) + 2 * BBOX_PAD_PX
        ax.add_patch(mpatches.Rectangle(
            (rect_x, rect_y), rect_w, rect_h,
            linewidth=1.5, edgecolor="white", facecolor="none",
        ))

        # White × at cell center
        ax.plot(cx, cy, marker="x", color="white", markersize=8, markeredgewidth=1.5)

        if row_i == 0:
            ax.set_title(COL_TITLES[col_i], fontsize=16, pad=6)

    axes[row_i][0].text(
        -0.12, 0.5, group,
        transform=axes[row_i][0].transAxes,
        fontsize=16, fontweight="bold",
        va="center", ha="right", rotation=90,
    )

plt.tight_layout(h_pad=0.3, w_pad=0.3)
plt.savefig(output_pdf, dpi=150, bbox_inches="tight")
print(f"Saved to {output_pdf}")
