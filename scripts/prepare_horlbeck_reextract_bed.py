#!/usr/bin/env python
# scripts/prepare_horlbeck_reextract_bed.py
"""
BED con las coordenadas de la matriz consolidada Horlbeck
para re-extraer H3K4me3 (y las demás) y comparar con los valores guardados.
"""

from pathlib import Path
import pandas as pd

MAT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
OUT_DIR = Path("data/interim/conditionA/reextract_horlbeck")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BED = OUT_DIR / "horlbeck_consolidated_coords.bed"
TABLE = OUT_DIR / "horlbeck_consolidated_for_reextract.csv"

def main():
    df = pd.read_csv(MAT)
    print(f"Matriz consolidada: {df.shape}")

    # Necesitamos chromosome + coordinate
    req = ["guide_sequence", "chromosome", "coordinate", "H3K4me3_mean",
           "ATAC_mean", "H3K27ac_mean", "H3K27me3_mean", "efficacy_score"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        print("Columnas disponibles:", df.columns.tolist())
        raise SystemExit(f"Faltan columnas: {missing}")

    df = df.dropna(subset=["chromosome", "coordinate"]).copy()
    df["coordinate"] = df["coordinate"].astype(int)
    df["start"] = df["coordinate"] - 1
    df["end"] = df["coordinate"]
    df["region_id"] = [f"h_{i:05d}" for i in range(len(df))]

    bed = df[["chromosome", "start", "end", "region_id"]].copy()
    bed["score"] = 0
    bed["strand"] = "."
    bed.to_csv(BED, sep="\t", header=False, index=False)

    df.to_csv(TABLE, index=False)
    print(f"Guías con coords válidas: {len(df):,}")
    print(f"BED:   {BED}")
    print(f"Tabla: {TABLE}")
    print(df[["guide_sequence", "chromosome", "coordinate", "H3K4me3_mean"]].head(3))

if __name__ == "__main__":
    main()