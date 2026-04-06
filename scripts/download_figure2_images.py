#!/usr/bin/env -S uv run --with pandas python3
"""Download Figure 2 representative images (all channels) from S3.

Usage: python scripts/download_figure2_images.py <csv_path> <output_dir>

Reads figure2_representative_cells.csv and downloads all three channels
(c1=DNA, c2=Actin, c3=Mito) for each representative image using a single
batched s5cmd run for maximum transfer concurrency.

S3 layout: .../Mito_Morphology_input/images/<subject> Mito_Morphology/<filename>
"""
import os
import subprocess
import sys
import tempfile

import pandas as pd

S3_BASE = "s3://imaging-platform/projects/2016_08_01_RadialMitochondriaDistribution_donna/Mito_Morphology_input/images"
CHANNELS = ["c1", "c2", "c3"]

csv_path, output_dir = sys.argv[1], sys.argv[2]
df = pd.read_csv(csv_path)

pairs = []
for _, row in df.iterrows():
    for channel in CHANNELS:
        filename = row["FileName_Mito"].replace("_c3_", f"_{channel}_")
        src = f"{S3_BASE}/{row['subject']} Mito_Morphology/{filename}"
        dst = f"{output_dir}/{filename}"
        pairs.append(f'cp "{src}" "{dst}"')

with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
    f.write("\n".join(pairs))
    batch_file = f.name

try:
    result = subprocess.run(["s5cmd", "run", batch_file], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    print(result.stdout.strip())
finally:
    os.unlink(batch_file)
