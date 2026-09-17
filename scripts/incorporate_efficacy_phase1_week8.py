#!/usr/bin/env python3

import os
import re
import glob
import numpy as np
import pandas as pd

BASE = "data/processed/sgRNA_feature_matrix_phase1.csv"
OUT = "data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv"
REPORT = "logs/week8_efficacy_report.md"

DATASET_DIRS = [
    "data/raw/crispr_datasets",
    "data/interim",
]

DATASETS = {
    "Horlbeck2016": ["Horlbeck2016", "horlbeck2016", "Horlbeck"],
    "sanson2018": ["Sanson2018", "sanson2018", "Sanson"],
    "Gasperini2019": ["Gasperini2019", "gasperini2019", "Gasperini"],
    "replogle2022": ["Replogle2022", "replogle2022", "Replogle"],
}

TARGET_PATTERNS = [
    "efficacy",
    "score_norm",
    "score",
    "log2fc",
    "log2_fc",
    "logfc",
    "phenotype",
    "activity",
    "effect",
    "fitness",
]

GUIDE_PATTERNS = [
    "sgrna_id",
    "sgRNA_id",
    "guide_id",
    "guide_sequence",
    "sgRNA",
    "guide",
]

print("=" * 70)
print("ChromaCRISPR Phase 1 — Week 8")
print("Paso 5 — Incorporación de variable de eficacia")
print("=" * 70)

# ------------------------------------------------------------------
# 1. Cargar matriz final
# ------------------------------------------------------------------

print("\nLoading Phase 1 matrix...")
base = pd.read_csv(BASE, low_memory=False)

print("Rows:", len(base))
print("Columns:", len(base.columns))

if "efficacy_score" in base.columns:
    raise RuntimeError(
        "La matriz ya contiene efficacy_score. "
        "No se debe sobrescribir sin inspección."
    )

# ------------------------------------------------------------------
# 2. Buscar archivos originales
# ------------------------------------------------------------------

print("\n=== DISCOVERING DATASETS ===")

files = []
for directory in DATASET_DIRS:
    if os.path.exists(directory):
        for root, _, filenames in os.walk(directory):
            for filename in filenames:
                path = os.path.join(root, filename)
                if filename.lower().endswith(
                    (".csv", ".tsv", ".txt", ".parquet")
                ):
                    files.append(path)

files = sorted(set(files))

print("Candidate tabular files:", len(files))

dataset_files = {name: [] for name in DATASETS}

for path in files:
    low = path.lower()

    for dataset, patterns in DATASETS.items():
        if any(pattern.lower() in low for pattern in patterns):
            dataset_files[dataset].append(path)

for dataset, paths in dataset_files.items():
    print(f"\n{dataset}:")
    for path in paths:
        print("  ", path)

# ------------------------------------------------------------------
# 3. Inspeccionar columnas sin asumir nombres
# ------------------------------------------------------------------

def read_table(path):
    ext = os.path.splitext(path)[1].lower()

    if ext == ".parquet":
        return pd.read_parquet(path)

    # Intentar autodetección para CSV/TSV/TXT
    try:
        return pd.read_csv(path, sep=None, engine="python", nrows=5)
    except Exception:
        return pd.read_csv(path, nrows=5)


def normalize_column(col):
    return re.sub(r"[^a-z0-9]+", "_", str(col).strip().lower()).strip("_")


def find_candidates(columns, patterns):
    result = []

    for col in columns:
        norm = normalize_column(col)

        for pattern in patterns:
            p = normalize_column(pattern)

            if norm == p or p in norm:
                result.append(col)
                break

    return result


metadata = {}

print("\n=== COLUMN INSPECTION ===")

for dataset, paths in dataset_files.items():

    metadata[dataset] = []

    for path in paths:
        try:
            sample = read_table(path)
        except Exception as e:
            print(f"\nWARNING: could not inspect {path}")
            print("Reason:", e)
            continue

        target_candidates = find_candidates(
            sample.columns,
            TARGET_PATTERNS
        )

        guide_candidates = find_candidates(
            sample.columns,
            GUIDE_PATTERNS
        )

        print(f"\n--- {dataset} ---")
        print("File:", path)
        print("Columns:")

        for i, col in enumerate(sample.columns):
            print(f"  {i}: {col}")

        print("Target candidates:", target_candidates)
        print("Guide/key candidates:", guide_candidates)

        metadata[dataset].append({
            "path": path,
            "columns": list(sample.columns),
            "target_candidates": target_candidates,
            "guide_candidates": guide_candidates,
        })

# ------------------------------------------------------------------
# 4. Detenerse si no hay candidatos suficientes
# ------------------------------------------------------------------

usable = {}

for dataset, entries in metadata.items():

    for entry in entries:

        if entry["target_candidates"] and entry["guide_candidates"]:
            usable.setdefault(dataset, []).append(entry)

print("\n=== USABLE DATASETS ===")

for dataset, entries in usable.items():
    print(dataset, ":", len(entries), "candidate file(s)")

if not usable:
    report = [
        "# ChromaCRISPR Phase 1 — Week 8",
        "",
        "## Paso 5 — Incorporación de eficacia",
        "",
        "No se encontró automáticamente ningún archivo con una combinación "
        "identificable de columna de guía y columna de eficacia/score.",
        "",
        "Se realizó inspección de los archivos candidatos encontrados en:",
        *[f"- `{x}`" for x in DATASET_DIRS],
        "",
        "El paso se detuvo sin modificar la matriz final.",
    ]

    with open(REPORT, "w") as fh:
        fh.write("\n".join(report) + "\n")

    print("\nPaso 5 detenido: no se identificó una llave y variable de eficacia.")
    print("Reporte:", REPORT)
    raise SystemExit(0)

# ------------------------------------------------------------------
# 5. Mostrar candidatos y NO asumir todavía una unión incorrecta
# ------------------------------------------------------------------

print("\n=== TARGET / KEY CANDIDATES ===")

for dataset, entries in usable.items():

    print(f"\n{dataset}")

    for entry in entries:
        print("File:", entry["path"])
        print("Target candidates:", entry["target_candidates"])
        print("Key candidates:", entry["guide_candidates"])

print(
    "\nIMPORTANTE: esta primera ejecución solamente identifica "
    "las columnas candidatas. No se realiza todavía una unión automática "
    "porque los nombres de score y las llaves deben verificarse contra "
    "los datasets originales."
)

# Guardar inventario para revisión
inventory = []

for dataset, entries in metadata.items():
    for entry in entries:
        inventory.append({
            "dataset": dataset,
            "file": entry["path"],
            "target_candidates": ";".join(
                map(str, entry["target_candidates"])
            ),
            "guide_candidates": ";".join(
                map(str, entry["guide_candidates"])
            ),
        })

inventory_df = pd.DataFrame(inventory)

report = [
    "# ChromaCRISPR Phase 1 — Week 8",
    "",
    "## Paso 5 — Identificación de variable de eficacia",
    "",
    f"- Matriz base: `{BASE}`",
    f"- Guías: {len(base):,}",
    "",
    "### Resultado",
    "",
    "Se inspeccionaron los archivos tabulares disponibles en los directorios "
    "de datos originales/intermedios.",
    "",
    "Se identificaron candidatos de columnas de eficacia y de llaves de "
    "unión, pero **no se realizó una integración automática en esta ejecución** "
    "para evitar asignar incorrectamente un score a una guía.",
    "",
    "### Candidatos detectados",
    "",
]

for _, row in inventory_df.iterrows():
    report.extend([
        f"#### {row['dataset']}",
        f"- Archivo: `{row['file']}`",
        f"- Targets candidatos: `{row['target_candidates']}`",
        f"- Llaves candidatas: `{row['guide_candidates']}`",
        "",
    ])

report.extend([
    "### Estado",
    "",
    "- Matriz `sgRNA_feature_matrix_phase1.csv`: no modificada.",
    "- `efficacy_score`: todavía no incorporado.",
    "- Se requiere verificar la correspondencia exacta entre cada dataset, "
      "su score y la llave `guide_sequence`/ID antes de efectuar el merge.",
])

with open(REPORT, "w") as fh:
    fh.write("\n".join(report) + "\n")

print("\nReport saved:", REPORT)
print("\nPaso 5 detenido después de la identificación de candidatos.")
