#!/bin/bash
set -euo pipefail

# ----------------------------
# User-configurable paths
# ----------------------------
FASTQ_DIR="./data/raw"
REF="./refs/hg38.fa"
PON=""
GERMLINE=""
INTERVALS=""

OUT_QC="./results/qc"
OUT_ALIGN="./results/alignment"
OUT_DEDUP="./results/umi_dedup"
OUT_MUTECT="./results/mutect"

THREADS=8

# ----------------------------
# 01 FASTQ QC
# ----------------------------
echo "[STEP 1] FASTQC"
python preprocessing/01_quality_control.py \
    --fastq_dir "$FASTQ_DIR" \
    --out_dir "$OUT_QC" \
    --threads "$THREADS"

# ----------------------------
# 02 Alignment
# ----------------------------
echo "[STEP 2] BWA-MEM2 alignment"
python preprocessing/02_alignment.py \
    --fastq_dir "$FASTQ_DIR" \
    --reference "$REF" \
    --out_dir "$OUT_ALIGN" \
    --threads "$THREADS"

# ----------------------------
# 03 UMI dedup
# ----------------------------
echo "[STEP 3] UMI-tools dedup"
python preprocessing/03_umi_processing_umitools.py \
    --bam_dir "$OUT_ALIGN" \
    --out_dir "$OUT_DEDUP" \
    --method directional \
    --threads "$THREADS"

# ----------------------------
# 04 Variant Calling (Mutect2)
# ----------------------------
echo "[STEP 4] Mutect2 variant calling"

CMD="python preprocessing/04_variant_calling_Mutect2.py \
    --bam_dir $OUT_DEDUP \
    --reference $REF \
    --out_dir $OUT_MUTECT \
    --threads $THREADS"

# Append optional parameters
if [[ -n "$PON" ]]; then
    CMD="$CMD --pon $PON"
fi

if [[ -n "$GERMLINE" ]]; then
    CMD="$CMD --germline $GERMLINE"
fi

if [[ -n "$INTERVALS" ]]; then
    CMD="$CMD --intervals $INTERVALS"
fi

echo "[INFO] Running: $CMD"
eval "$CMD"

echo "[COMPLETE] ctDNA pipeline finished successfully."
