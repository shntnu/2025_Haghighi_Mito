#!/usr/bin/env -S uv run --with tifffile --with matplotlib --with numpy --with pandas python3
"""Compose Figure 2: representative single cells per diagnosis group.

Usage: python scripts/compose_figure2.py <cells_csv> <images_dir> <output_pdf>

Layout: 5 rows (BP, Control, MDD, SZ, SZA) × 5 columns (DNA, Mito, Actin, Merge, Zoomed).
Canonical subjects matching upstream: MCL128, MCL162, 287, MCL015, 263.

Channels: c1=DNA (blue), c2=Actin (green), c3=Mito (red).
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
GROUP_ORDER = ["BP", "Control", "MDD", "SZ", "SZA"]

# Channel display: c1=DNA/blue, c2=Actin/green, c3=Mito/red
CHANNEL_COLOR = {"c1": np.array([0.2, 0.4, 1.0]), "c2": np.array([0.2, 0.9, 0.2]), "c3": np.array([1.0, 0.2, 0.2])}
CHANNEL_LABEL = {"c1": "DNA", "c2": "Actin", "c3": "Mito"}
COL_ORDER = ["c1", "c3", "c2", "merge", "zoom"]
COL_TITLES = ["DNA", "Mito", "Actin", "Merge", "Cell (zoom)"]

ZOOM_PAD_PX = 30  # pixels of padding around bounding box for zoom crop


def normalize(img, plo=1, phi=99):
    """Clip to percentiles and scale to [0, 1]."""
    lo, hi = np.percentile(img, plo), np.percentile(img, phi)
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1)


def make_merge(channels: dict) -> np.ndarray:
    """Combine normalized single-channel images into an RGB composite."""
    h, w = next(iter(channels.values())).shape
    rgb = np.zeros((h, w, 3), dtype=float)
    for ch, img_norm in channels.items():
        rgb += img_norm[:, :, None] * CHANNEL_COLOR[ch]
    return np.clip(rgb, 0, 1)


def load_channels(images_dir: Path, filename_mito: str) -> dict:
    """Load all three channels for one image, returning normalized arrays."""
    result = {}
    for ch in ("c1", "c2", "c3"):
        fname = filename_mito.replace("_c3_", f"_{ch}_")
        path = images_dir / fname
        raw = tifffile.imread(str(path))
        result[ch] = normalize(raw)
    return result


def crop(img, rmin, rmax, cmin, cmax, pad=ZOOM_PAD_PX):
    """Crop image with optional padding, clamped to image bounds."""
    h, w = img.shape[:2]
    r0 = max(0, rmin - pad)
    r1 = min(h, rmax + pad)
    c0 = max(0, cmin - pad)
    c1 = min(w, cmax + pad)
    return img[r0:r1, c0:c1]


cells_csv, images_dir, output_pdf = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
df = pd.read_csv(cells_csv)
df["subject"] = df["subject"].astype(str)

fig, axes = plt.subplots(
    len(GROUP_ORDER), len(COL_ORDER),
    figsize=(5 * len(COL_ORDER), 4 * len(GROUP_ORDER)),
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

    # Bounding box (BoundingBoxMinimum/Maximum are x=col, y=row in CellProfiler)
    bx_min, by_min = int(r["Cells_AreaShape_BoundingBoxMinimum_X"]), int(r["Cells_AreaShape_BoundingBoxMinimum_Y"])
    bx_max, by_max = int(r["Cells_AreaShape_BoundingBoxMaximum_X"]), int(r["Cells_AreaShape_BoundingBoxMaximum_Y"])

    for col_i, col_key in enumerate(COL_ORDER):
        ax = axes[row_i][col_i]
        ax.axis("off")

        if col_key == "merge":
            ax.imshow(merge, interpolation="nearest")
            # Draw box around selected cell
            rect = mpatches.Rectangle(
                (bx_min, by_min), bx_max - bx_min, by_max - by_min,
                linewidth=1.5, edgecolor="white", facecolor="none",
            )
            ax.add_patch(rect)
        elif col_key == "zoom":
            zoomed = crop(merge, by_min, by_max, bx_min, bx_max)
            ax.imshow(zoomed, interpolation="nearest")
        else:
            ax.imshow(channels[col_key], cmap="gray", interpolation="nearest")

        if row_i == 0:
            ax.set_title(COL_TITLES[col_i], fontsize=11, pad=4)

    # Row label: add text annotation to the left of the first column
    axes[row_i][0].text(
        -0.08, 0.5, group,
        transform=axes[row_i][0].transAxes,
        fontsize=12, fontweight="bold",
        va="center", ha="right", rotation=90,
    )

plt.tight_layout(h_pad=0.5, w_pad=0.5)
plt.savefig(output_pdf, dpi=150, bbox_inches="tight")
print(f"Saved to {output_pdf}")
