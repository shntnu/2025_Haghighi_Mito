#!/usr/bin/env -S uv run --with tifffile --with matplotlib --with numpy --with pandas python3
"""Compose Figure 2: representative single cells per diagnosis group.

Usage: python scripts/compose_figure2.py <cells_csv> <images_dir> <output_pdf>

Layout: 5 rows (Control, BP, SZ, SZA, MDD) × 6 columns
        (DNA, Mito, Actin, composite, Mito crop color, Mito crop grayscale).
Canonical subjects matching upstream: MCL162, MCL128, MCL015, 263, 287.

Channels: c1=DNA (blue), c2=Actin (green), c3=Mito (red).
Normalization: 0.5/99.5 percentile for full-field columns;
               min/max for grayscale zoom (matches figure legend).
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
COL_ORDER = ["c1", "c3", "c2", "merge", "zoom_merge", "zoom_gray"]
COL_TITLES = ["DNA", "Mito", "Actin", "composite", "", ""]

ZOOM_PAD_PX = 30


def normalize(img, plo=0.5, phi=99.5):
    lo, hi = np.percentile(img, plo), np.percentile(img, phi)
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return np.clip((img.astype(float) - lo) / (hi - lo), 0, 1)


def normalize_minmax(img):
    lo, hi = float(img.min()), float(img.max())
    if hi == lo:
        return np.zeros_like(img, dtype=float)
    return (img.astype(float) - lo) / (hi - lo)


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


def crop_box(img, rmin, rmax, cmin, cmax, pad=ZOOM_PAD_PX):
    h, w = img.shape[:2]
    return img[max(0, rmin - pad):min(h, rmax + pad),
               max(0, cmin - pad):min(w, cmax + pad)]


def pad_to_height(img, target_h):
    """Embed img centered vertically in a black canvas of target_h rows."""
    h = img.shape[0]
    if h >= target_h:
        return img
    pad_top = (target_h - h) // 2
    pad_bot = target_h - h - pad_top
    pad_width = ((pad_top, pad_bot), (0, 0)) if img.ndim == 2 else ((pad_top, pad_bot), (0, 0), (0, 0))
    return np.pad(img, pad_width)


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

    bx_min = int(r["Cells_AreaShape_BoundingBoxMinimum_X"])
    by_min = int(r["Cells_AreaShape_BoundingBoxMinimum_Y"])
    bx_max = int(r["Cells_AreaShape_BoundingBoxMaximum_X"])
    by_max = int(r["Cells_AreaShape_BoundingBoxMaximum_Y"])

    full_h = channels["c3"].shape[0]
    mito_crop_color = pad_to_height(crop_box(merge, by_min, by_max, bx_min, bx_max), full_h)
    mito_crop_gray = pad_to_height(normalize_minmax(crop_box(channels["c3"], by_min, by_max, bx_min, bx_max)), full_h)

    for col_i, col_key in enumerate(COL_ORDER):
        ax = axes[row_i][col_i]
        ax.axis("off")

        if col_key == "merge":
            ax.imshow(merge, interpolation="nearest")
            ax.add_patch(mpatches.Rectangle(
                (bx_min, by_min), bx_max - bx_min, by_max - by_min,
                linewidth=1.5, edgecolor="white", facecolor="none",
            ))
        elif col_key == "zoom_merge":
            ax.imshow(mito_crop_color, interpolation="nearest")
        elif col_key == "zoom_gray":
            ax.imshow(mito_crop_gray, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
        else:
            ax.imshow(channels[col_key], cmap="gray", interpolation="nearest")

        if row_i == 0:
            ax.set_title(COL_TITLES[col_i], fontsize=11, pad=4)

    axes[row_i][0].text(
        -0.08, 0.5, group,
        transform=axes[row_i][0].transAxes,
        fontsize=12, fontweight="bold",
        va="center", ha="right", rotation=90,
    )

plt.tight_layout(h_pad=0.3, w_pad=0.3)
plt.savefig(output_pdf, dpi=150, bbox_inches="tight")
print(f"Saved to {output_pdf}")
