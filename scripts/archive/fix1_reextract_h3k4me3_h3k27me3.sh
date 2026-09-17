#!/bin/bash
# scripts/fix1_reextract_h3k4me3_h3k27me3.sh
# Re-extrae H3K4me3 y H3K27me3 sobre sgRNA_hg38.bed
# con los mismos parámetros que ATAC/H3K27ac (±500 bp, bin 10).

set -euo pipefail

BED="data/interim/features/sgRNA_hg38.bed"
OUTDIR="data/interim/features/histones_fixed"
mkdir -p "$OUTDIR"
BW_DIR="data/raw/encode_k562"
WINDOW=500
BINSIZE=10

echo "=== Verificando BED ==="
wc -l "$BED"
head -3 "$BED"

for MARK in H3K4me3 H3K27me3; do
    case $MARK in
        H3K4me3)  BW="${BW_DIR}/H3K4me3_K562.bigWig" ;;
        H3K27me3) BW="${BW_DIR}/H3K27me3_K562.bigWig" ;;
    esac

    echo ""
    echo "=== $MARK ==="
    computeMatrix reference-point \
        --scoreFileName "$BW" \
        --regionsFileName "$BED" \
        --referencePoint center \
        --beforeRegionStartLength $WINDOW \
        --afterRegionStartLength $WINDOW \
        --binSize $BINSIZE \
        --missingDataAsZero \
        --outFileName "${OUTDIR}/${MARK}_matrix_fixed.gz" \
        --outFileNameMatrix "${OUTDIR}/${MARK}_matrix_fixed.tab" \
        --numberOfProcessors 6 \
        2> "${OUTDIR}/${MARK}_computeMatrix.log"

    echo "$MARK terminado."
    ls -lah "${OUTDIR}/${MARK}_matrix_fixed.gz"
done

echo ""
echo "=== LISTO Fix-1 (matrices) ==="
ls -lah "$OUTDIR"