#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

MATRIX = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
XLSX = Path("data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx")

matrix = pd.read_csv(MATRIX, low_memory=False)
crispri = pd.read_excel(XLSX, sheet_name="CRISPRi")

# Normalización
matrix["guide_sequence"] = (
    matrix["guide_sequence"].astype("string").str.strip().str.upper()
)
matrix["gene_symbol"] = (
    matrix["gene_symbol"].astype("string").str.strip().str.upper()
)

crispri["sgRNA sequence"] = (
    crispri["sgRNA sequence"].astype("string").str.strip().str.upper()
)
crispri["gene symbol"] = (
    crispri["gene symbol"].astype("string").str.strip().str.upper()
)

key = ["guide_sequence", "gene_symbol"]

activity_keys = (
    crispri[["sgRNA sequence", "gene symbol"]]
    .drop_duplicates()
    .rename(
        columns={
            "sgRNA sequence": "guide_sequence",
            "gene symbol": "gene_symbol",
        }
    )
)

# Match de todas las filas de la matriz contra las llaves Horlbeck
matched = matrix.merge(
    activity_keys,
    on=key,
    how="inner"
)

print("=" * 80)
print("DIAGNÓSTICO DEL SUBSET HORLBECK2016")
print("=" * 80)

print(f"\nFilas matriz principal              : {len(matrix):,}")
print(f"Llaves únicas Horlbeck             : {len(activity_keys):,}")
print(f"Filas matriz que hacen match       : {len(matched):,}")
print(
    f"Filas matriz sin match             : "
    f"{len(matrix) - len(matched):,}"
)

print("\nDistribución de datasets entre los matches:")
print(
    matched["dataset"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nDistribución dataset + experiment:")
print(
    matched.groupby(
        ["dataset", "experiment"],
        dropna=False
    ).size().to_string()
)

# ¿Cuántos matches son exactamente filas Horlbeck?
horlbeck_rows = matched[
    matched["dataset"].astype("string").str.lower() == "horlbeck2016"
]

print(f"\nFilas dataset=Horlbeck2016          : {len(horlbeck_rows):,}")

# Comparar número de filas Horlbeck con source
print(f"Filas source CRISPRi                : {len(crispri):,}")

# Duplicación de llaves dentro de Horlbeck
activity_dup = (
    crispri.groupby(
        ["sgRNA sequence", "gene symbol"],
        dropna=False
    )
    .size()
    .reset_index(name="n_rows")
)

activity_dup = activity_dup[activity_dup["n_rows"] > 1]

print(
    f"\nLlaves duplicadas dentro de Horlbeck: "
    f"{len(activity_dup):,}"
)

if len(activity_dup):
    print("\nPRIMERAS DUPLICADAS EN HORLBECK:")
    print(activity_dup.head(20).to_string(index=False))

# Verificar consistencia de scores para cada llave
score_col = "CRISPRi activity score [Horlbeck et al., eLife 2016]"

score_check = (
    crispri.groupby(
        ["sgRNA sequence", "gene symbol"],
        dropna=False
    )[score_col]
    .agg(["count", "nunique"])
    .reset_index()
)

inconsistent = score_check[score_check["nunique"] > 1]

print(
    f"\nLlaves con scores diferentes para la misma "
    f"secuencia+gen: {len(inconsistent):,}"
)

if len(inconsistent):
    print("\nPRIMERAS INCONSISTENCIAS:")
    print(inconsistent.head(20).to_string(index=False))

# Ver si todas las filas Horlbeck de la matriz tienen llave en source
matrix_h = matrix[
    matrix["dataset"].astype("string").str.lower() == "horlbeck2016"
]

matrix_h_keys = matrix_h[key].drop_duplicates()

missing_h = matrix_h_keys.merge(
    activity_keys,
    on=key,
    how="left",
    indicator=True
)

missing_h = missing_h[
    missing_h["_merge"] == "left_only"
]

print(
    f"\nLlaves Horlbeck en matriz sin source match: "
    f"{len(missing_h):,}"
)

print("\n" + "=" * 80)
print("DIAGNÓSTICO COMPLETADO")
print("=" * 80)
