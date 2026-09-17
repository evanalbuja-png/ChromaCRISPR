#!/bin/bash
# scripts/reextract_horlbeck_chromatin.sh

set -euo pipefail

BED="data/interim/conditionA/reextract_horlbeck/horlbeck_consolidated_coords.bed"
OUTDIR="data/interim/conditionA/reextract_horlbeck/matrices"
mkdir -p "$OUTDIR"
BW_DIR="data/raw/encode_k562"
WINDOW=500
BINSIZE=10

for MARK in H3K4me3 ATAC H3K27ac H3K27me3; do
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
        --outFileName "${OUTDIR}/${MARK}_reextract.gz" \
        --outFileNameMatrix "${OUTDIR}/${MARK}_reextract.tab" \
        --numberOfProcessors 6 \
        2> "${OUTDIR}/${MARK}_reextract.log"
    echo "$MARK terminado."
done

ls -lah "$OUTDIR"
echo "=== LISTO R2 ==="