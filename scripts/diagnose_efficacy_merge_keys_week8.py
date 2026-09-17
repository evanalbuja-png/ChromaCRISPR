#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

MATRIX = Path("data/processed/sgRNA_feature_matrix_phase1.csv")

print("=" * 80)
print("DIAGNÓSTICO DE LLAVES PARA INTEGRACIÓN CRISPRi EFFICACY")
print("=" * 80)

matrix = pd.read_csv(MATRIX, low_memory=False)

matrix["guide_sequence"] = (
    matrix["guide_sequence"].astype("string").str.strip().str.upper()
)

matrix["gene_symbol"] = (
    matrix["gene_symbol"].astype("string").str.strip().str.upper()
)

key = ["guide_sequence", "gene_symbol"]

# Duplicados de la llave
dup = (
    matrix.groupby(key, dropna=False)
    .size()
    .reset_index(name="n_rows")
)

dup = dup[dup["n_rows"] > 1].sort_values(
    "n_rows", ascending=False
)

print(f"\nTotal de filas matriz: {len(matrix):,}")
print(f"Llaves únicas: {matrix[key].drop_duplicates().shape[0]:,}")
print(f"Llaves duplicadas: {len(dup):,}")
print(f"Filas pertenecientes a llaves duplicadas: {dup['n_rows'].sum():,}")

print("\nDistribución del número de filas por llave duplicada:")
print(
    dup["n_rows"]
    .value_counts()
    .sort_index()
    .to_string()
)

print("\nPrimeras 30 llaves duplicadas:")
print(dup.head(30).to_string(index=False))

# Mostrar contexto completo de las primeras llaves duplicadas
if len(dup) > 0:
    print("\n" + "=" * 80)
    print("CONTEXTO DE LAS PRIMERAS 10 LLAVES DUPLICADAS")
    print("=" * 80)

    for _, row in dup.head(10).iterrows():
        seq = row["guide_sequence"]
        gene = row["gene_symbol"]

        subset = matrix[
            (matrix["guide_sequence"] == seq) &
            (matrix["gene_symbol"] == gene)
        ]

        print("\n---")
        print(f"guide_sequence : {seq}")
        print(f"gene_symbol    : {gene}")
        print(subset[
            [
                "sgrna_id",
                "dataset",
                "experiment",
                "guide_sequence",
                "gene_symbol",
                "gene_id",
                "chromosome",
                "coordinate",
                "strand"
            ]
        ].to_string(index=False))

print("\n" + "=" * 80)
print("DIAGNÓSTICO COMPLETADO")
print("=" * 80)
