#!/bin/bash
###############################################################################
# Full QIIME2 pipeline for the temperature synthetic community 16S rRNA dataset.
#
# Workflow:
#   1. Cutadapt primer removal (both batches, paired-end)
#   2. Re-summarize quality after primer removal
#   3. DADA2 paired-end denoising (both batches, separately)
#   4. QC visualizations
#   5. Merge tables and representative sequences across batches
#   6. Export feature table and ASV sequences
#
# Region    : V4-V5 (515F-907R)
# Primers   : 515F GTGCCAGCMGCCGCGGTAA  (19 nt)
#             907R CCGTCAATTCMTTTRAGTTT (20 nt)
# Trunc     : 230 (R1) / 205 (R2) -- consistent with prior R-DADA2 pipeline
#
# Taxonomy classification is intentionally not included; ASVs will be
# annotated by local BLAST against the known synthetic-community isolates.
###############################################################################

set -eo pipefail

# ============================================================================
# Paths and parameters
# ============================================================================
WORK_DIR="${WORK_DIR:-/lustre/home/zfhu/16s_project/analysis}"
cd "${WORK_DIR}"

PRIMER_F="GTGCCAGCMGCCGCGGTAA"   # 515F
PRIMER_R="CCGTCAATTCMTTTRAGTTT"  # 907R

TRUNC_LEN_F=228
TRUNC_LEN_R=205
MAX_EE_F=2
MAX_EE_R=2
N_THREADS=${N_THREADS:-32}

# ============================================================================
# Activate QIIME2 environment
# ============================================================================
source "${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh}"
conda activate "${QIIME2_ENV:-qiime2-amplicon-2023.9}"

echo "===================================================================="
echo "QIIME2 version : $(qiime --version 2>&1 | head -1)"
echo "Start time     : $(date)"
echo "Working dir    : ${WORK_DIR}"
echo "Threads        : ${N_THREADS}"
echo "Trunc lengths  : F=${TRUNC_LEN_F}, R=${TRUNC_LEN_R}"
echo "===================================================================="

# ============================================================================
# Pre-flight checks
# ============================================================================
for f in demux_old.qza demux_new.qza metadata.tsv; do
    if [ ! -f "${WORK_DIR}/${f}" ]; then
        echo "ERROR: ${f} not found in ${WORK_DIR}"
        exit 1
    fi
done

# ============================================================================
# Step 1: Primer removal with cutadapt
# ============================================================================
for batch in old new; do
    if [ ! -f "${WORK_DIR}/demux_${batch}_trimmed.qza" ]; then
        echo "[$(date '+%H:%M:%S')] Removing primers from batch ${batch}..."
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences "${WORK_DIR}/demux_${batch}.qza" \
            --p-front-f "${PRIMER_F}" \
            --p-front-r "${PRIMER_R}" \
            --p-discard-untrimmed \
            --p-cores ${N_THREADS} \
            --o-trimmed-sequences "${WORK_DIR}/demux_${batch}_trimmed.qza" \
            --verbose
    else
        echo "[$(date '+%H:%M:%S')] demux_${batch}_trimmed.qza exists, skipping cutadapt."
    fi
done

# ============================================================================
# Step 2: Summarize quality after primer removal
# ============================================================================
echo "[$(date '+%H:%M:%S')] Summarizing trimmed sequence quality..."
for batch in old new; do
    qiime demux summarize \
        --i-data "${WORK_DIR}/demux_${batch}_trimmed.qza" \
        --o-visualization "${WORK_DIR}/demux_${batch}_trimmed.qzv"
done

# ============================================================================
# Step 3: DADA2 denoising (paired-end), per batch
# ============================================================================
for batch in old new; do
    if [ ! -f "${WORK_DIR}/table_${batch}.qza" ]; then
        echo "[$(date '+%H:%M:%S')] Denoising batch ${batch} with DADA2..."
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "${WORK_DIR}/demux_${batch}_trimmed.qza" \
            --p-trim-left-f  0 \
            --p-trim-left-r  0 \
            --p-trunc-len-f  ${TRUNC_LEN_F} \
            --p-trunc-len-r  ${TRUNC_LEN_R} \
            --p-max-ee-f     ${MAX_EE_F} \
            --p-max-ee-r     ${MAX_EE_R} \
            --p-n-threads    ${N_THREADS} \
            --o-table                    "${WORK_DIR}/table_${batch}.qza" \
            --o-representative-sequences "${WORK_DIR}/rep-seqs_${batch}.qza" \
            --o-denoising-stats          "${WORK_DIR}/denoising-stats_${batch}.qza" \
            --verbose
    else
        echo "[$(date '+%H:%M:%S')] table_${batch}.qza exists, skipping DADA2."
    fi
done

# ============================================================================
# Step 4: QC visualizations
# ============================================================================
echo "[$(date '+%H:%M:%S')] Generating QC visualizations..."
for batch in old new; do
    qiime metadata tabulate \
        --m-input-file "${WORK_DIR}/denoising-stats_${batch}.qza" \
        --o-visualization "${WORK_DIR}/denoising-stats_${batch}.qzv"
done

# ============================================================================
# Step 5: Merge tables and representative sequences across batches
# ============================================================================
if [ ! -f "${WORK_DIR}/table_merged.qza" ]; then
    echo "[$(date '+%H:%M:%S')] Merging feature tables across batches..."
    qiime feature-table merge \
        --i-tables "${WORK_DIR}/table_old.qza" \
        --i-tables "${WORK_DIR}/table_new.qza" \
        --o-merged-table "${WORK_DIR}/table_merged.qza"

    qiime feature-table merge-seqs \
        --i-data "${WORK_DIR}/rep-seqs_old.qza" \
        --i-data "${WORK_DIR}/rep-seqs_new.qza" \
        --o-merged-data "${WORK_DIR}/rep-seqs_merged.qza"
fi

# ============================================================================
# Step 6: Summarize merged outputs
# ============================================================================
echo "[$(date '+%H:%M:%S')] Summarizing merged outputs..."
qiime feature-table summarize \
    --i-table "${WORK_DIR}/table_merged.qza" \
    --m-sample-metadata-file "${WORK_DIR}/metadata.tsv" \
    --o-visualization "${WORK_DIR}/table_merged.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${WORK_DIR}/rep-seqs_merged.qza" \
    --o-visualization "${WORK_DIR}/rep-seqs_merged.qzv"

# ============================================================================
# Step 7: Export final outputs
# ============================================================================
echo "[$(date '+%H:%M:%S')] Exporting final tables..."
mkdir -p "${WORK_DIR}/exported"

qiime tools export \
    --input-path  "${WORK_DIR}/table_merged.qza" \
    --output-path "${WORK_DIR}/exported"

biom convert \
    -i "${WORK_DIR}/exported/feature-table.biom" \
    -o "${WORK_DIR}/exported/feature-table.tsv" \
    --to-tsv

qiime tools export \
    --input-path  "${WORK_DIR}/rep-seqs_merged.qza" \
    --output-path "${WORK_DIR}/exported"

echo "===================================================================="
echo "End time : $(date)"
echo "All steps completed successfully."
echo "===================================================================="
echo ""
echo "Final outputs in ${WORK_DIR}/exported/:"
echo "  feature-table.tsv     - ASV abundance table"
echo "  dna-sequences.fasta   - ASV representative sequences"
echo ""
echo "Taxonomy: annotate ASVs by local BLAST against the known"
echo "isolate 16S sequences (manual step, not part of this pipeline)."
echo ""
echo "Inspect visualizations on https://view.qiime2.org :"
echo "  demux_old_trimmed.qzv / demux_new_trimmed.qzv  - post-primer-trim QC"
echo "  denoising-stats_old.qzv / denoising-stats_new.qzv  - per-sample retention"
echo "  table_merged.qzv      - merged ASV table summary"
echo "  rep-seqs_merged.qzv   - ASV sequence summary"
