#!/usr/bin/env python
# scripts/diagnose_coordinate_offset.py
"""
Diagnosticar el origen del offset de coordenadas entre Horlbeck y Sanson
para las mismas guide_sequence.
"""

from pathlib import Path
import pandas as pd
import numpy as np

HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
SANSON = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
SANSON_REF = Path("data/processed/sanson2018_reference_FINAL.csv")
# Si tienes la referencia Horlbeck original:
HORLBECK_REF_CANDIDATES = [
    "data/processed/horlbeck2016_reference_FINAL.csv",
    "data/processed/sgRNA_unified_FINAL.csv",
    "data/processed/sgRNA_feature_reference_hg38.csv",
]

def main():
    h = pd.read_csv(HORLBECK)
    s = pd.read_csv(SANSON)
    m = h.merge(s, on="guide_sequence", how="inner", suffixes=("_h", "_s"))
    print(f"Overlap: {len(m):,}\n")

    # Diferencia de coordenadas
    m["delta"] = m["coordinate_h"] - m["coordinate_s"]
    print("=== Distribución del delta (coord_h - coord_s) ===")
    print(m["delta"].describe())
    print(f"\nValores más frecuentes de delta:")
    print(m["delta"].value_counts().head(15))

    # ¿Depende del strand?
    if "strand_h" in m.columns or "strand_s" in m.columns:
        strand_col = "strand_s" if "strand_s" in m.columns else "strand_h"
        print(f"\n=== Delta por strand ({strand_col}) ===")
        print(m.groupby(strand_col)["delta"].agg(["mean", "median", "std", "count"]))

    # Mirar la referencia Sanson original
    print("\n=== Referencia Sanson (cómo se asignó coordinate) ===")
    ref = pd.read_csv(SANSON_REF)
    print(ref.columns.tolist())
    print(ref[["guide_sequence", "chromosome", "coordinate", "strand", "offset", "status"]].head(3))

    # Intentar cargar referencia Horlbeck
    print("\n=== Buscando referencia Horlbeck ===")
    for cand in HORLBECK_REF_CANDIDATES:
        p = Path(cand)
        if p.exists():
            print(f"Encontrado: {cand}")
            hr = pd.read_csv(p, nrows=5)
            print("Columnas:", hr.columns.tolist())
            # buscar columnas de coordenada
            for c in hr.columns:
                if any(x in c.lower() for x in ["coord", "start", "pos", "chrom"]):
                    print(f"  {c}: {hr[c].head(2).tolist()}")
            break
    else:
        print("No se encontró referencia Horlbeck clara.")

    # Muestra de 5 guías con delta pequeño y 5 con delta grande
    print("\n=== Ejemplos delta ≈ 0–30 bp ===")
    close = m[m["delta"].abs() <= 30].head(5)
    cols = ["guide_sequence", "coordinate_h", "coordinate_s", "delta"]
    for c in ["strand_h", "strand_s", "chromosome_h", "chromosome_s"]:
        if c in m.columns:
            cols.append(c)
    print(close[cols].to_string())

    print("\n=== Ejemplos |delta| > 1000 bp ===")
    far = m[m["delta"].abs() > 1000].head(5)
    print(far[cols].to_string() if len(far) else "Ninguno")

if __name__ == "__main__":
    main()