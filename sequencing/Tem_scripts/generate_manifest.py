#!/usr/bin/env python3
"""
Generate QIIME2 manifest files and sample metadata for the temperature
synthetic community 16S rRNA dataset.

Experimental design
-------------------
- Synthetic communities: 12 communities under 5 temperature conditions (W1-W5)
- Time course: day 1-10
- Replicates:
    * R1: full time course (day 1-10)
    * R2 / R3: day 8, day 9, day 10 only

Sequencing batches
------------------
- Batch 1 (Temperature_V1): plates P1-P8
    * R1 time course (day 1-10) for all communities under W1-W5
    * R2 / R3 endpoint (day 10) for all communities under W1-W5
- Batch 2 (Temperature_V2): plates P1-P3
    * R2 / R3 for day 8 (P2) and day 9 (P1) under W1-W4
    * R2 / R3 for day 8 and day 9 under W5 (P3)

Output
------
- manifest_old.tsv : QIIME2 PairedEndFastqManifestPhred33V2 for batch 1
- manifest_new.tsv : QIIME2 PairedEndFastqManifestPhred33V2 for batch 2
- metadata.tsv     : sample metadata for all samples
"""

import os
from collections import Counter

# ============================================================================
# Paths
# ============================================================================
OLD_RAW_DIR = os.environ.get("OLD_RAW_DIR", "/lustre/home/zfhu/16s_project/Temperature_V1")
NEW_RAW_DIR = os.environ.get("NEW_RAW_DIR", "/lustre/home/zfhu/16s_project/Temperature_V2")
OUT_DIR = os.environ.get("OUT_DIR", "/lustre/home/zfhu/16s_project/analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================
# Batch 1 well-to-sample mapping (plates P1-P8)
# ============================================================================
# Plates P1-P8 contain community samples (P9 isolate plate is excluded from
# this analysis). Samples are arranged sequentially: each set of 12 wells
# corresponds to one community, and every 12 communities cycle through one
# temperature condition in the order W1, W2, W3, W4, W5.
#
# Within each 12-well community block:
#   columns 1-10 : day 1-10, replicate R1
#   column   11  : day 10,   replicate R2
#   column   12  : day 10,   replicate R3
# ----------------------------------------------------------------------------
CONDITIONS_OLD = ["W1", "W2", "W3", "W4", "W5"]

def decode_old_community(plate, num):
    """Map (plate, well_number) on P1-P8 to (day, condition, community, replica)."""
    global_num    = (plate - 1) * 96 + num
    community_idx = (global_num - 1) // 12
    col           = (global_num - 1) %  12 + 1
    condition_idx = community_idx // 12
    condition     = CONDITIONS_OLD[condition_idx]
    community     = community_idx % 12 + 1
    if col <= 10:
        day, replica = col, 1
    elif col == 11:
        day, replica = 10, 2
    else:
        day, replica = 10, 3
    return day, condition, community, replica

# ============================================================================
# Batch 2 well-to-sample mapping (plates P1-P3)
# ============================================================================
# P1 holds day 9 samples; P2 holds day 8 samples. Each plate is divided into
# four 6-column x 4-row quadrants:
#     top-left    (col 1-6,  row 1-4)  -> W1
#     top-right   (col 7-12, row 1-4)  -> W3
#     bottom-left (col 1-6,  row 5-8)  -> W2
#     bottom-right(col 7-12, row 5-8)  -> W4
#
# Within each quadrant (local_col 1-6, local_row 1-4):
#   community = 2 * local_col - 1  if local_row in {1, 2}
#             = 2 * local_col      if local_row in {3, 4}
#   replica   = 2  if local_row in {1, 3}
#             = 3  if local_row in {2, 4}
#
# P3 supplies the W5 condition for day 8 (col 1) and day 9 (col 2), rows 1-6.
# ----------------------------------------------------------------------------
PLATE_TO_DAY_NEW  = {"P1": 9, "P2": 8}
PART_TO_CONDITION = {
    (0, 0): "W1",
    (0, 4): "W2",
    (6, 0): "W3",
    (6, 4): "W4",
}
P3_COL_TO_DAY = {1: 8, 2: 9}
P3_ROW_MAP    = {
    1: (3,  2),  2: (3,  3),
    3: (7,  2),  4: (7,  3),
    5: (12, 2),  6: (12, 3),
}

def decode_new_p12(plate_str, col, row):
    """Map (plate, col, row) on P1/P2 to (day, condition, community, replica)."""
    day        = PLATE_TO_DAY_NEW[plate_str]
    col_offset = 6 if col > 6 else 0
    row_offset = 4 if row > 4 else 0
    condition  = PART_TO_CONDITION.get((col_offset, row_offset))
    if condition is None:
        return None
    local_col = col - col_offset
    local_row = row - row_offset
    community = 2 * local_col - 1 if local_row in (1, 2) else 2 * local_col
    replica   = 2 if local_row in (1, 3) else 3
    return day, condition, community, replica

def decode_new_p3(col, row):
    """Map (col, row) on P3 to (day, condition, community, replica) for W5."""
    day      = P3_COL_TO_DAY.get(col)
    comm_rep = P3_ROW_MAP.get(row)
    if day is None or comm_rep is None:
        return None
    community, replica = comm_rep
    return day, "W5", community, replica

# ============================================================================
# FASTQ file discovery
# ============================================================================
def find_fastq_old(plate, num):
    """Locate paired-end FASTQ files for a batch 1 sample directory."""
    sample_dir = os.path.join(OLD_RAW_DIR, f"P{plate}-{num}")
    if not os.path.exists(sample_dir):
        return None, None
    files = os.listdir(sample_dir)
    r1 = sorted(f for f in files if "_R1_" in f and f.endswith(".fastq.gz"))
    r2 = sorted(f for f in files if "_R2_" in f and f.endswith(".fastq.gz"))
    if not r1 or not r2:
        return None, None
    return os.path.join(sample_dir, r1[0]), os.path.join(sample_dir, r2[0])

def find_fastq_new(orig_name):
    """Locate paired-end FASTQ files for a batch 2 sample directory."""
    sample_dir = os.path.join(NEW_RAW_DIR, orig_name)
    if not os.path.exists(sample_dir):
        return None, None
    files = os.listdir(sample_dir)
    r1 = sorted(f for f in files if ".R1." in f and f.endswith(".fq.gz"))
    r2 = sorted(f for f in files if ".R2." in f and f.endswith(".fq.gz"))
    if not r1 or not r2:
        return None, None
    return os.path.join(sample_dir, r1[0]), os.path.join(sample_dir, r2[0])

# ============================================================================
# Build sample records
# ============================================================================
records = []
missing = []

# ---- Batch 1: community samples on plates P1-P8 ----------------------------
print("Building batch 1 records (plates P1-P8)...")
for plate in range(1, 9):
    max_num = 48 if plate == 8 else 96
    for num in range(1, max_num + 1):
        result = decode_old_community(plate, num)
        if result is None:
            continue
        day, condition, community, replica = result
        r1, r2 = find_fastq_old(plate, num)
        if r1 is None:
            missing.append(f"P{plate}-{num}")
            continue
        records.append({
            "sample_id": f"temperature_{condition}_C{community:02d}_R{replica}_D{day}",
            "r1": r1, "r2": r2,
            "condition": condition, "community": community,
            "replica": replica, "day": day,
            "batch": "old", "orig": f"P{plate}-{num}",
        })

# ---- Batch 2: P1 and P2 ----------------------------------------------------
print("Building batch 2 records (plates P1, P2)...")
for plate_str in ["P1", "P2"]:
    for col in range(1, 13):
        for row in range(1, 9):
            result = decode_new_p12(plate_str, col, row)
            if result is None:
                continue
            day, condition, community, replica = result
            orig = f"{plate_str}-{col}-{row}"
            r1, r2 = find_fastq_new(orig)
            if r1 is None:
                missing.append(orig)
                continue
            records.append({
                "sample_id": f"temperature_{condition}_C{community:02d}_R{replica}_D{day}",
                "r1": r1, "r2": r2,
                "condition": condition, "community": community,
                "replica": replica, "day": day,
                "batch": "new", "orig": orig,
            })

# ---- Batch 2: P3 (W5 supplement) -------------------------------------------
print("Building batch 2 records (plate P3, W5 supplement)...")
for col in range(1, 3):
    for row in range(1, 7):
        result = decode_new_p3(col, row)
        if result is None:
            continue
        day, condition, community, replica = result
        orig = f"P3-{col}-{row}"
        r1, r2 = find_fastq_new(orig)
        if r1 is None:
            missing.append(orig)
            continue
        records.append({
            "sample_id": f"temperature_{condition}_C{community:02d}_R{replica}_D{day}",
            "r1": r1, "r2": r2,
            "condition": condition, "community": community,
            "replica": replica, "day": day,
            "batch": "new", "orig": orig,
        })

# ============================================================================
# Write outputs
# ============================================================================
records.sort(key=lambda x: (x["condition"], x["community"], x["replica"], x["day"]))

old_records = [r for r in records if r["batch"] == "old"]
new_records = [r for r in records if r["batch"] == "new"]

for batch_name, batch_records in [("old", old_records), ("new", new_records)]:
    out_path = os.path.join(OUT_DIR, f"manifest_{batch_name}.tsv")
    with open(out_path, "w") as f:
        f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for r in batch_records:
            f.write(f"{r['sample_id']}\t{r['r1']}\t{r['r2']}\n")

with open(os.path.join(OUT_DIR, "metadata.tsv"), "w") as f:
    f.write("sample-id\tcondition\tcommunity\treplica\tday\tbatch\n")
    for r in records:
        f.write(f"{r['sample_id']}\t{r['condition']}\t{r['community']}\t"
                f"{r['replica']}\t{r['day']}\t{r['batch']}\n")

# ============================================================================
# Summary
# ============================================================================
print(f"\nTotal samples written : {len(records)}")
print(f"  Batch 1 (old)       : {len(old_records)}")
print(f"  Batch 2 (new)       : {len(new_records)}")
print(f"Missing FASTQ files   : {len(missing)}")
if missing:
    for m in missing[:10]:
        print(f"    {m}")
    if len(missing) > 10:
        print(f"    ... and {len(missing) - 10} more")

cond_counts = Counter(r["condition"] for r in records)
print(f"By condition          : {dict(sorted(cond_counts.items()))}")