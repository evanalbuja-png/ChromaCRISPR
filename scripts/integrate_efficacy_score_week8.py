#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

MATRIX = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
XLSX = Path("data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx")
OUTPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv")
LOG = Path("logs/week8_log.md")

# 1. Cargar datos
matrix = pd.read_csv(MATRIX)
crispri = pd.read_excel(XLSX, sheet_name="CRISPRi")

# Normalizar nombres/secuencias para la llave
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

score_col = "CRISPRi activity score [Horlbeck et al., eLife 2016]"

# 2. Preparar tabla de efficacy usando la llave secuencia + gen
efficacy = crispri[
    ["sgRNA sequence", "gene symbol", score_col]
].copy()

efficacy = efficacy.rename(
    columns={
        "sgRNA sequence": "guide_sequence",
        "gene symbol": "gene_symbol",
        score_col: "efficacy_score",
    }
)

# Verificar que no existan múltiples scores para una misma llave
duplicates = efficacy.duplicated(
    subset=["guide_sequence", "gene_symbol"],
    keep=False
)

if duplicates.any():
    dup_count = duplicates.sum()
    raise ValueError(
        f"ERROR: {dup_count} registros CRISPRi pertenecen a "
        "llaves guide_sequence + gene_symbol duplicadas."
    )

# 3. Unión left: conserva las 155,454 guías de la matriz principal
merged = matrix.merge(
    efficacy,
    on=["guide_sequence", "gene_symbol"],
    how="left",
    validate="one_to_one",
)

# 4. Validaciones
if len(merged) != len(matrix):
    raise ValueError(
        f"ERROR: la integración cambió el número de filas: "
        f"{len(matrix)} -> {len(merged)}"
    )

matched = merged["efficacy_score"].notna().sum()
total = len(merged)
coverage = matched / total * 100

# 5. Estadísticas
score_stats = merged["efficacy_score"].dropna()

print("=" * 80)
print("INTEGRACIÓN CRISPRi ACTIVITY SCORE")
print("=" * 80)
print(f"Matriz original             : {len(matrix):,}")
print(f"Matriz integrada            : {len(merged):,}")
print(f"Guías con efficacy_score    : {matched:,}")
print(f"Guías sin efficacy_score    : {total - matched:,}")
print(f"Cobertura                   : {coverage:.2f}%")
print()
print("Distribución efficacy_score")
print(f"Mean                        : {score_stats.mean():.12g}")
print(f"Median                      : {score_stats.median():.12g}")
print(f"Min                         : {score_stats.min():.12g}")
print(f"Max                         : {score_stats.max():.12g}")
print("=" * 80)

# 6. Guardar matriz actualizada
merged.to_csv(OUTPUT, index=False)

# 7. Actualizar log
LOG.parent.mkdir(parents=True, exist_ok=True)

with LOG.open("a", encoding="utf-8") as f:
    f.write("\n## Integration of CRISPRi efficacy score\n\n")
    f.write(
        "- Input matrix: "
        "`data/processed/sgRNA_feature_matrix_phase1.csv`\n"
    )
    f.write(
        "- CRISPRi source: "
        "`data/raw/crispr_datasets/Horlbeck2016/"
        "elife-19760-supp1-v2.xlsx` (sheet `CRISPRi`)\n"
    )
    f.write("- Join key: `guide_sequence + gene_symbol`\n")
    f.write(f"- Total guides: {total:,}\n")
    f.write(f"- Guides with efficacy_score: {matched:,}\n")
    f.write(f"- Guides without efficacy_score: {total - matched:,}\n")
    f.write(f"- Coverage: {coverage:.2f}%\n")
    f.write(f"- Mean: {score_stats.mean():.12g}\n")
    f.write(f"- Median: {score_stats.median():.12g}\n")
    f.write(f"- Min: {score_stats.min():.12g}\n")
    f.write(f"- Max: {score_stats.max():.12g}\n")
    f.write(
        "- Output: "
        "`data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv`\n"
    )

print(f"\nOutput guardado en: {OUTPUT}")
print(f"Log actualizado: {LOG}")
