#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

MATRIX = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
XLSX = Path("data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx")
OUTPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv")
LOG = Path("logs/week8_log.md")

KEY = ["guide_sequence", "gene_symbol"]
SCORE_COL = "CRISPRi activity score [Horlbeck et al., eLife 2016]"

# 1. Cargar datos
matrix = pd.read_csv(MATRIX, low_memory=False)
activity = pd.read_excel(XLSX, sheet_name="CRISPRi")

# El Excel usa nombres diferentes
activity = activity.rename(columns={
    "sgRNA sequence": "guide_sequence",
    "gene symbol": "gene_symbol"
})

# Verificar columnas requeridas
required_activity = KEY + [SCORE_COL]
missing = [c for c in required_activity if c not in activity.columns]

if missing:
    raise ValueError(
        f"ERROR: faltan columnas en Horlbeck CRISPRi: {missing}"
    )

# 2. Normalizar llaves
matrix["guide_sequence"] = (
    matrix["guide_sequence"]
    .astype("string")
    .str.strip()
    .str.upper()
)

matrix["gene_symbol"] = (
    matrix["gene_symbol"]
    .astype("string")
    .str.strip()
)

activity["guide_sequence"] = (
    activity["guide_sequence"]
    .astype("string")
    .str.strip()
    .str.upper()
)

activity["gene_symbol"] = (
    activity["gene_symbol"]
    .astype("string")
    .str.strip()
)

# 3. Verificar unicidad del archivo Horlbeck
activity_key_dups = activity.duplicated(KEY).sum()

if activity_key_dups != 0:
    raise ValueError(
        f"ERROR: Horlbeck CRISPRi contiene "
        f"{activity_key_dups} llaves duplicadas."
    )

# 4. Inicializar columna
matrix["efficacy_score"] = pd.NA

# 5. Seleccionar únicamente Horlbeck2016
horlbeck_mask = matrix["dataset"].eq("Horlbeck2016")
horlbeck = matrix.loc[horlbeck_mask, KEY].copy()

# 6. Merge controlado
merged = horlbeck.merge(
    activity[KEY + [SCORE_COL]],
    on=KEY,
    how="left",
    validate="many_to_one"
)

# 7. Asignar score solamente a Horlbeck2016
matrix.loc[horlbeck_mask, "efficacy_score"] = (
    merged[SCORE_COL].to_numpy()
)

# 8. Métricas
n_total = len(matrix)
n_horlbeck = horlbeck_mask.sum()
n_scored = matrix["efficacy_score"].notna().sum()

coverage = (
    100 * n_scored / n_horlbeck
    if n_horlbeck > 0 else 0
)

scores = pd.to_numeric(
    matrix["efficacy_score"],
    errors="coerce"
).dropna()

print("=" * 80)
print("INTEGRACIÓN CONTROLADA DEL CRISPRi ACTIVITY SCORE")
print("=" * 80)

print(f"Filas matriz total            : {n_total:,}")
print(f"Filas Horlbeck2016            : {n_horlbeck:,}")
print(f"Guías con efficacy_score      : {n_scored:,}")
print(f"Cobertura Horlbeck2016        : {coverage:.2f}%")
print()

print("Distribución del efficacy_score:")
print(f"Mean   : {scores.mean():.12f}")
print(f"Median : {scores.median():.12f}")
print(f"Min    : {scores.min():.12f}")
print(f"Max    : {scores.max():.12f}")
print()

# 9. Comprobaciones
non_horlbeck_scored = (
    (~horlbeck_mask) &
    matrix["efficacy_score"].notna()
).sum()

if non_horlbeck_scored != 0:
    raise ValueError(
        f"ERROR: {non_horlbeck_scored} filas fuera de Horlbeck2016 "
        "recibieron efficacy_score."
    )

# 10. Guardar matriz
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
matrix.to_csv(OUTPUT, index=False)

print(f"Output guardado: {OUTPUT}")

# 11. Actualizar log
LOG.parent.mkdir(parents=True, exist_ok=True)

with LOG.open("a", encoding="utf-8") as f:
    f.write("\n## Integración controlada CRISPRi efficacy score\n")
    f.write("- Fuente: Horlbeck2016, hoja CRISPRi\n")
    f.write("- Llave: guide_sequence + gene_symbol\n")
    f.write("- Restricción: dataset == Horlbeck2016\n")
    f.write(f"- Filas matriz: {n_total:,}\n")
    f.write(f"- Filas Horlbeck2016: {n_horlbeck:,}\n")
    f.write(f"- efficacy_score no nulo: {n_scored:,}\n")
    f.write(f"- Cobertura: {coverage:.2f}%\n")
    f.write(f"- Mean: {scores.mean():.12f}\n")
    f.write(f"- Median: {scores.median():.12f}\n")
    f.write(f"- Min: {scores.min():.12f}\n")
    f.write(f"- Max: {scores.max():.12f}\n")
    f.write(f"- Filas no-Horlbeck con score: {non_horlbeck_scored}\n")
    f.write(f"- Output: {OUTPUT}\n")
