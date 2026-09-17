#!/bin/bash
# scripts/extract_chromatin_overlap_horlbeck_coords.sh

set -euo pipefail

BED="data/interim/conditionA/overlap_control/overlap_horlbeck_coords.bed"
OUTDIR="data/interim/conditionA/overlap_control/matrices"
mkdir -p "$OUTDIR"
BW_DIR="data/raw/encode_k562"
WINDOW=500
BINSIZE=10

for MARK in ATAC H3K27ac H3K4me3 H3K27me3; do
    case $MARK in
        ATAC)     BW="${BW_DIR}/ATAC-seq_K562.bigWig" ;;
        H3K27ac)  BW="${BW_DIR}/H3K27ac_K562.bigWig" ;;
        H3K4me3)  BW="${BW_DIR}/H3K4me3_K562.bigWig" ;;
        H3K27me3) BW="${BW_DIR}/H3K27me3_K562.bigWig" ;;
    esac

    echo "=== $MARK ==="
    computeMatrix reference-point \
        --scoreFileName "$BW" \
        --regionsFileName "$BED" \
        --referencePoint center \
        --beforeRegionStartLength $WINDOW \
        --afterRegionStartLength $WINDOW \
        --binSize $BINSIZE \
        --skipZeros \
        --outFileName "${OUTDIR}/${MARK}_overlap_matrix.gz" \
        --outFileNameMatrix "${OUTDIR}/${MARK}_overlap_matrix.tab" \
        --numberOfProcessors 4 \
        2> "${OUTDIR}/${MARK}_overlap.log"
    echo "$MARK terminado."
done

echo ""
ls -lah "$OUTDIR"
echo "=== LISTO C2 ==="