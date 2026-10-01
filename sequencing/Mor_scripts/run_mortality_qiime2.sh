#!/bin/bash
###############################################################################
# Full QIIME2 pipeline for the mortality 16S rRNA dataset.
#
# Workflow:
#   1. Import single-end primer-trimmed merged FASTQ files
#   2. Summarize sequence quality
#   3. DADA2 denoising (single-end)
#   4. Export QC visualizations
#   5. Taxonomic classification with SILVA 138 V4 pre-trained classifier
#   6. Export feature table, representative sequences, and taxonomy
#
# Region : V4 (515F-806R), primers already removed, reads already merged
# Trunc  : 250 bp (based on read length distribution; primary mode at 252 bp)
###############################################################################

set -eo pipefail

# ============================================================================
# Paths and parameters
# ============================================================================
WORK_DIR="${WORK_DIR:-/lustre/home/zfhu/Mortality/analysis}"
REF_DIR="${REF_DIR:-/lustre/home/zfhu/Mortality/reference}"
CLASSIFIER="${REF_DIR}/silva-138-99-515-806-nb-classifier.qza"
cd "${WORK_DIR}"

TRIM_LEFT=0
TRUNC_LEN=250
MAX_EE=2
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
echo "Classifier     : ${CLASSIFIER}"
echo "Threads        : ${N_THREADS}"
echo "Trunc length   : ${TRUNC_LEN}"
echo "===================================================================="

# ============================================================================
# Pre-flight checks
# ============================================================================
if [ ! -f "${WORK_DIR}/manifest_communities.tsv" ]; then
    echo "ERROR: manifest_communities.tsv not found in ${WORK_DIR}"
    exit 1
fi

if [ ! -f "${CLASSIFIER}" ]; then
    echo "ERROR: classifier not found at ${CLASSIFIER}"
    exit 1
fi

# ============================================================================
# Step 1: Import sequences
# ============================================================================
if [ ! -f "${WORK_DIR}/demux.qza" ]; then
    echo "[$(date '+%H:%M:%S')] Importing sequences..."
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path  "${WORK_DIR}/manifest_communities.tsv" \
        --output-path "${WORK_DIR}/demux.qza" \
        --input-format SingleEndFastqManifestPhred33V2
else
    echo "[$(date '+%H:%M:%S')] demux.qza exists, skipping import."
fi

# ============================================================================
# Step 2: Summarize sequence quality
# ============================================================================
echo "[$(date '+%H:%M:%S')] Summarizing demultiplexed sequences..."
qiime demux summarize \
    --i-data "${WORK_DIR}/demux.qza" \
    --o-visualization "${WORK_DIR}/demux.qzv"

# ============================================================================
# Step 3: DADA2 denoising
# ============================================================================
if [ ! -f "${WORK_DIR}/table.qza" ]; then
    echo "[$(date '+%H:%M:%S')] Denoising with DADA2 (trunc-len=${TRUNC_LEN})..."
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs "${WORK_DIR}/demux.qza" \
        --p-trim-left  ${TRIM_LEFT} \
        --p-trunc-len  ${TRUNC_LEN} \
        --p-max-ee     ${MAX_EE} \
        --p-n-threads  ${N_THREADS} \
        --o-table                    "${WORK_DIR}/table.qza" \
        --o-representative-sequences "${WORK_DIR}/rep-seqs.qza" \
        --o-denoising-stats          "${WORK_DIR}/denoising-stats.qza" \
        --verbose
else
    echo "[$(date '+%H:%M:%S')] table.qza exists, skipping denoising."
fi

# ============================================================================
# Step 4: QC visualizations
# ============================================================================
echo "[$(date '+%H:%M:%S')] Generating QC visualizations..."

qiime metadata tabulate \
    --m-input-file "${WORK_DIR}/denoising-stats.qza" \
    --o-visualization "${WORK_DIR}/denoising-stats.qzv"

qiime feature-table summarize \
    --i-table "${WORK_DIR}/table.qza" \
    --m-sample-metadata-file "${WORK_DIR}/metadata.tsv" \
    --o-visualization "${WORK_DIR}/table.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${WORK_DIR}/rep-seqs.qza" \
    --o-visualization "${WORK_DIR}/rep-seqs.qzv"

# ============================================================================
# Step 5: Taxonomic classification
# ============================================================================
if [ ! -f "${WORK_DIR}/taxonomy.qza" ]; then
    echo "[$(date '+%H:%M:%S')] Classifying ASVs against SILVA 138 (V4)..."
    qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER}" \
        --i-reads "${WORK_DIR}/rep-seqs.qza" \
        --p-n-jobs ${N_THREADS} \
        --o-classification "${WORK_DIR}/taxonomy.qza"
else
    echo "[$(date '+%H:%M:%S')] taxonomy.qza exists, skipping classification."
fi

qiime metadata tabulate \
    --m-input-file "${WORK_DIR}/taxonomy.qza" \
    --o-visualization "${WORK_DIR}/taxonomy.qzv"

# ============================================================================
# Step 6: Export final results
# ============================================================================
echo "[$(date '+%H:%M:%S')] Exporting final tables..."
mkdir -p "${WORK_DIR}/exported"

qiime tools export \
    --input-path  "${WORK_DIR}/table.qza" \
    --output-path "${WORK_DIR}/exported"

biom convert \
    -i "${WORK_DIR}/exported/feature-table.biom" \
    -o "${WORK_DIR}/exported/feature-table.tsv" \
    --to-tsv

qiime tools export \
    --input-path  "${WORK_DIR}/rep-seqs.qza" \
    --output-path "${WORK_DIR}/exported"

qiime tools export \
    --input-path  "${WORK_DIR}/taxonomy.qza" \
    --output-path "${WORK_DIR}/exported"

echo "===================================================================="
echo "End time : $(date)"
echo "All steps completed successfully."
echo "===================================================================="
echo ""
echo "Final outputs in ${WORK_DIR}/exported/:"
echo "  feature-table.tsv     - ASV abundance table"
echo "  dna-sequences.fasta   - ASV representative sequences"
echo "  taxonomy.tsv          - ASV taxonomic assignments"
echo ""
echo "Inspect visualizations on https://view.qiime2.org :"
echo "  demux.qzv             - input sequence quality summary"
echo "  denoising-stats.qzv   - per-sample retention rates"
echo "  table.qzv             - feature table summary"
echo "  rep-seqs.qzv          - ASV sequence summary"
echo "  taxonomy.qzv          - classification results"
