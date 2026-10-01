#!/usr/bin/env python3
"""
Generate QIIME2 manifest and metadata for the mortality 16S rRNA dataset.

Experimental design
-------------------
- Synthetic communities: 12 communities under 5 conditions (W1-W5)
- Time course: day 1-6, replicate R1
- 16S region: V4 (515F-806R), reads are pre-merged and primer-trimmed

Input filename convention
-------------------------
"W-C-D.clean.fastq.gz", e.g. "5-12-6.clean.fastq.gz" denotes
condition W5, community 12, day 6.

Output
------
- manifest_communities.tsv : QIIME2 SingleEndFastqManifestPhred33V2
- metadata.tsv             : sample metadata
"""

import os
from collections import Counter

# ============================================================================
# Paths
# ============================================================================
BASE_DIR = os.environ.get("BASE_DIR", "/lustre/home/zfhu/Mortality")
COMM_DIR = os.path.join(BASE_DIR, "Communities")
OUT_DIR  = os.path.join(BASE_DIR, "analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================
# Build sample records
# ============================================================================
print("Scanning community samples...")
records = []

for fname in sorted(os.listdir(COMM_DIR)):
    if not fname.endswith(".clean.fastq.gz"):
        continue

    name = fname.replace(".clean.fastq.gz", "")
    w, c, d = map(int, name.split("-"))

    sample_id = f"mortality_W{w}_C{c:02d}_R1_D{d}"
    abs_path  = os.path.join(COMM_DIR, fname)
    records.append({
        "sample_id": sample_id,
        "filepath":  abs_path,
        "condition": f"W{w}",
        "community": c,
        "replica":   1,
        "day":       d,
    })

# ============================================================================
# Sort and write outputs
# ============================================================================
records.sort(key=lambda x: (x["condition"], x["community"], x["day"]))

manifest_path = os.path.join(OUT_DIR, "manifest_communities.tsv")
with open(manifest_path, "w") as f:
    f.write("sample-id\tabsolute-filepath\n")
    for r in records:
        f.write(f"{r['sample_id']}\t{r['filepath']}\n")

metadata_path = os.path.join(OUT_DIR, "metadata.tsv")
with open(metadata_path, "w") as f:
    f.write("sample-id\tcondition\tcommunity\treplica\tday\n")
    for r in records:
        f.write(f"{r['sample_id']}\t{r['condition']}\t{r['community']}\t"
                f"{r['replica']}\t{r['day']}\n")

# ============================================================================
# Summary
# ============================================================================
print(f"\nTotal samples written : {len(records)}")

cc = Counter(r["condition"] for r in records)
print(f"\nBy condition:")
for c in sorted(cc):
    print(f"  {c}: {cc[c]}")

dc = Counter(r["day"] for r in records)
print(f"\nBy day:")
for d in sorted(dc):
    print(f"  D{d}: {dc[d]}")

print(f"\nOutput files:")
print(f"  {manifest_path}")
print(f"  {metadata_path}")
