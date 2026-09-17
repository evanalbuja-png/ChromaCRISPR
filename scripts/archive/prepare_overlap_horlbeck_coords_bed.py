#!/usr/bin/env python
# scripts/prepare_overlap_horlbeck_coords_bed.py
"""
Genera BED de las 526 guías compartidas usando las coordenadas de Horlbeck
(para re-extraer features de cromatina en el anclaje del entrenamiento).
"""

from pathlib import Path
import pandas as pd

HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
SANSON = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
OUT_DIR = Path("data/interim/conditionA/overlap_control")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BED_OUT = OUT_DIR / "overlap_horlbeck_coords.bed"
TABLE_OUT = OUT_DIR / "overlap_guides_horlbeck_coords.csv"

def main():
    h = pd.read_csv(HORLBECK)
    s = pd.read_csv(SANSON)

    # Solo columnas necesarias de Horlbeck
    h_sub = h[["guide_sequence", "chromosome", "coordinate", "efficacy_score"]].copy()
    h_sub = h_sub.rename(columns={
        "chromosome": "chromosome_h",
        "coordinate": "coordinate_h"
    })

    # Traer LFC de Sanson
    s_sub = s[["guide_sequence", "lfc_HT29", "lfc_A375", "gene_symbol"]].copy()

    m = h_sub.merge(s_sub, on="guide_sequence", how="inner")
    print(f"Overlap con efficacy_score + LFC: {len(m):,}")

    # Limpiar
    m = m.dropna(subset=["chromosome_h", "coordinate_h"]).copy()
    m["coordinate_h"] = m["coordinate_h"].astype(int)
    m["start"] = m["coordinate_h"] - 1
    m["end"] = m["coordinate_h"]
    m["region_id"] = [f"ov_{i:04d}" for i in range(len(m))]
    m["strand"] = "."   # no lo necesitamos para reference-point center

    # BED
    bed = m[["chromosome_h", "start", "end", "region_id"]].copy()
    bed.columns = ["chrom", "start", "end", "name"]
    bed["score"] = 0
    bed["strand"] = "."
    bed.to_csv(BED_OUT, sep="\t", header=False, index=False)

    m.to_csv(TABLE_OUT, index=False)
    print(f"BED:   {BED_OUT}  ({len(bed)} regiones)")
    print(f"Tabla: {TABLE_OUT}")
    print(m[["guide_sequence", "chromosome_h", "coordinate_h", "efficacy_score", "lfc_HT29"]].head(3))

if __name__ == "__main__":
    main()