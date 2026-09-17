#!/usr/bin/env python
# scripts/prepare_sanson_bed_conditionA.py
"""
Prepara las guías Sanson2018 con coordenadas hg38 válidas
para extracción de features K562 (Condition A).
"""

from pathlib import Path
import pandas as pd
import numpy as np

REF = Path("data/processed/sanson2018_reference_FINAL.csv")
LFC = Path("data/processed/sanson2018_lfc_scores.csv")
OUT_DIR = Path("data/interim/conditionA")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BED_OUT = OUT_DIR / "sanson_hg38.bed"
TABLE_OUT = OUT_DIR / "sanson_guides_with_coords.csv"

def main():
    print("=== Cargando Sanson reference + LFC ===")
    ref = pd.read_csv(REF, low_memory=False)
    lfc = pd.read_csv(LFC, low_memory=False)

    print(f"Reference shape: {ref.shape}")
    print(f"LFC shape:       {lfc.shape}")

    # Unir por secuencia (las dos tablas tienen 114061 filas y están alineadas por orden)
    # Verificamos que las secuencias coincidan
    assert (ref["guide_sequence"] == lfc["sgRNA_sequence"]).all(), "Secuencias no alineadas!"

    df = ref.copy()
    df["lfc_HT29"] = lfc["lfc_HT29"].values
    df["lfc_A375"] = lfc["lfc_A375"].values
    df["gene_id"] = lfc["gene_id"].values
    df["set"] = lfc["set"].values

    # Filtrar solo guías con coordenadas válidas
    mask = df["chromosome"].notna() & df["coordinate"].notna()
    df_valid = df.loc[mask].copy()
    print(f"\nGuías con coordenadas válidas: {len(df_valid):,} / {len(df):,}")
    print(f"Guías sin coordenadas (descartadas): {(~mask).sum():,}")

    # Construir BED (0-based half-open para deepTools / bedtools)
    # coordinate parece ser el centro o el start del protospacer.
    # En el pipeline original se usó una ventana centrada → tratamos coordinate como centro.
    df_valid["coordinate"] = df_valid["coordinate"].astype(int)
    df_valid["start"] = df_valid["coordinate"] - 1          # 0-based start (aprox)
    df_valid["end"]   = df_valid["coordinate"]              # 1-based end → half-open

    # Orden canónico de cromosomas
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    df_valid["chromosome"] = pd.Categorical(df_valid["chromosome"], categories=chrom_order, ordered=True)
    df_valid = df_valid.sort_values(["chromosome", "start"]).reset_index(drop=True)

    # region_id único
    df_valid["region_id"] = [f"sanson_{i:06d}" for i in range(len(df_valid))]

    # Escribir BED (chrom, start, end, name, score, strand)
    bed = df_valid[["chromosome", "start", "end", "region_id"]].copy()
    bed["score"] = 0
    bed["strand"] = df_valid["strand"].fillna(".").values
    bed.to_csv(BED_OUT, sep="\t", header=False, index=False)

    # Tabla completa
    cols_keep = [
        "region_id", "guide_sequence", "gene_symbol", "gene_id", "set",
        "chromosome", "coordinate", "start", "end", "strand",
        "coordinate_source", "status", "offset",
        "lfc_HT29", "lfc_A375"
    ]
    df_valid[cols_keep].to_csv(TABLE_OUT, index=False)

    print(f"\nBED escrito:   {BED_OUT}  ({len(bed):,} regiones)")
    print(f"Tabla escrita: {TABLE_OUT}")
    print("\nDistribución por cromosoma (top 10):")
    print(df_valid["chromosome"].value_counts().head(10))
    print("\n=== LISTO Paso 4 ===")

if __name__ == "__main__":
    main()