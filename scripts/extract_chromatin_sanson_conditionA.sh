#!/bin/bash
# scripts/extract_chromatin_sanson_conditionA.sh
# Extrae matrices ATAC + H3K27ac + H3K4me3 + H3K27me3 de K562
# sobre las 100.942 coordenadas de Sanson (Condition A)

set -euo pipefail

BED="data/interim/conditionA/sanson_hg38.bed"
OUTDIR="data/interim/conditionA/matrices"
mkdir -p "$OUTDIR"

BW_DIR="data/raw/encode_k562"

# Ventana idéntica a la usada en el corpus de entrenamiento
WINDOW=500          # ±500 bp
BINSIZE=10

echo "=== 1. ATAC-seq K562 ==="
computeMatrix reference-point \
    --scoreFileName "${BW_DIR}/ATAC-seq_K562.bigWig" \
    --regionsFileName "$BED" \
    --referencePoint center \
    --beforeRegionStartLength $WINDOW \
    --afterRegionStartLength $WINDOW \
    --binSize $BINSIZE \
    --skipZeros \
    --outFileName "${OUTDIR}/atac_sanson_matrix.gz" \
    --outFileNameMatrix "${OUTDIR}/atac_sanson_matrix.tab" \
    --numberOfProcessors 8 \
    2> "${OUTDIR}/atac_sanson_computeMatrix.log"

echo "ATAC terminado."

echo "=== 2. H3K27ac K562 ==="
computeMatrix reference-point \
    --scoreFileName "${BW_DIR}/H3K27ac_K562.bigWig" \
    --regionsFileName "$BED" \
    --referencePoint center \
    --beforeRegionStartLength $WINDOW \
    --afterRegionStartLength $WINDOW \
    --binSize $BINSIZE \
    --skipZeros \
    --outFileName "${OUTDIR}/h3k27ac_sanson_matrix.gz" \
    --outFileNameMatrix "${OUTDIR}/h3k27ac_sanson_matrix.tab" \
    --numberOfProcessors 8 \
    2> "${OUTDIR}/h3k27ac_sanson_computeMatrix.log"

echo "H3K27ac terminado."

echo "=== 3. H3K4me3 K562 ==="
computeMatrix reference-point \
    --scoreFileName "${BW_DIR}/H3K4me3_K562.bigWig" \
    --regionsFileName "$BED" \
    --referencePoint center \
    --beforeRegionStartLength $WINDOW \
    --afterRegionStartLength $WINDOW \
    --binSize $BINSIZE \
    --skipZeros \
    --outFileName "${OUTDIR}/h3k4me3_sanson_matrix.gz" \
    --outFileNameMatrix "${OUTDIR}/h3k4me3_sanson_matrix.tab" \
    --numberOfProcessors 8 \
    2> "${OUTDIR}/h3k4me3_sanson_computeMatrix.log"

echo "H3K4me3 terminado."

echo "=== 4. H3K27me3 K562 ==="
computeMatrix reference-point \
    --scoreFileName "${BW_DIR}/H3K27me3_K562.bigWig" \
    --regionsFileName "$BED" \
    --referencePoint center \
    --beforeRegionStartLength $WINDOW \
    --afterRegionStartLength $WINDOW \
    --binSize $BINSIZE \
    --skipZeros \
    --outFileName "${OUTDIR}/h3k27me3_sanson_matrix.gz" \
    --outFileNameMatrix "${OUTDIR}/h3k27me3_sanson_matrix.tab" \
    --numberOfProcessors 8 \
    2> "${OUTDIR}/h3k27me3_sanson_computeMatrix.log"

echo "H3K27me3 terminado."

echo ""
echo "=== Resumen de archivos generados ==="
ls -lah "$OUTDIR"
echo ""
echo "=== LISTO Paso 5 ==="
echo "Revisa los .log si hay muchos 'Skipping regions'."