#!/usr/bin/env python3

import os
import re
import pandas as pd

SEARCH_DIRS = [
    "data/raw/crispr_datasets",
    "data/interim",
]

DATASETS = {
    "Horlbeck2016": [
        "horlbeck2016",
        "horlbeck",
    ],
    "Gasperini2019": [
        "gasperini2019",
        "gasperini",
    ],
    "Replogle2022": [
        "replogle2022",
        "replogle",
    ],
}

# Palabras que pueden identificar variables de eficacia/actividad.
SCORE_TERMS = [
    "score",
    "efficacy",
    "efficiency",
    "activity",
    "phenotype",
    "fitness",
    "log2fc",
    "log2_fc",
    "logfc",
    "fold",
    "effect",
    "activity",
    "depletion",
    "enrichment",
]

# Posibles llaves.
KEY_TERMS = [
    "sgrna",
    "sg_rna",
    "guide",
    "sequence",
    "grna",
    "id",
]

REPORT = "logs/week8_original_efficacy_inspection.md"


def normalize(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")


def read_table(path):
    ext = os.path.splitext(path)[1].lower()

    if ext == ".parquet":
        return pd.read_parquet(path)

    if ext in [".csv"]:
        return pd.read_csv(path, low_memory=False, nrows=5)

    if ext in [".tsv", ".txt"]:
        try:
            return pd.read_csv(
                path,
                sep="\t",
                low_memory=False,
                nrows=5,
            )
        except Exception:
            return pd.read_csv(
                path,
                sep=None,
                engine="python",
                low_memory=False,
                nrows=5,
            )

    return None


def find_columns(columns, terms):
    hits = []

    for col in columns:
        n = normalize(col)

        for term in terms:
            if normalize(term) in n:
                hits.append(col)
                break

    return hits


print("=" * 70)
print("ChromaCRISPR Phase 1 — Week 8")
print("Paso 5A — Inspección de datasets originales de eficacia")
print("=" * 70)

all_files = []

for directory in SEARCH_DIRS:
    if not os.path.exists(directory):
        continue

    for root, _, filenames in os.walk(directory):
        for filename in filenames:

            if filename.startswith("."):
                continue

            ext = os.path.splitext(filename)[1].lower()

            if ext in [".csv", ".tsv", ".txt", ".parquet"]:
                all_files.append(os.path.join(root, filename))

all_files = sorted(set(all_files))

print("\nTotal candidate files:", len(all_files))

results = []

for dataset, patterns in DATASETS.items():

    print("\n" + "=" * 70)
    print(dataset)
    print("=" * 70)

    matched_files = []

    for path in all_files:
        low = path.lower()

        if any(pattern in low for pattern in patterns):
            matched_files.append(path)

    if not matched_files:
        print("NO FILE MATCHED")
        results.append({
            "dataset": dataset,
            "file": "NONE",
            "score_candidates": "",
            "key_candidates": "",
        })
        continue

    for path in matched_files:

        print("\n--- FILE ---")
        print(path)

        try:
            df = read_table(path)
        except Exception as e:
            print("READ ERROR:", e)

            results.append({
                "dataset": dataset,
                "file": path,
                "score_candidates": "READ_ERROR",
                "key_candidates": "",
            })

            continue

        if df is None:
            continue

        print("\nColumns:")

        for i, col in enumerate(df.columns):
            print(f"{i}: {col}")

        score_candidates = find_columns(
            df.columns,
            SCORE_TERMS
        )

        key_candidates = find_columns(
            df.columns,
            KEY_TERMS
        )

        print("\nPotential SCORE/EFFICACY columns:")
        if score_candidates:
            for col in score_candidates:
                print("  -", col)
        else:
            print("  NONE")

        print("\nPotential GUIDE/KEY columns:")
        if key_candidates:
            for col in key_candidates:
                print("  -", col)
        else:
            print("  NONE")

        results.append({
            "dataset": dataset,
            "file": path,
            "score_candidates": "; ".join(
                map(str, score_candidates)
            ),
            "key_candidates": "; ".join(
                map(str, key_candidates)
            ),
        })

# ------------------------------------------------------------------
# Save report
# ------------------------------------------------------------------

report = []

report.append("# ChromaCRISPR Phase 1 — Week 8")
report.append("")
report.append("## Paso 5A — Inspección de datasets originales de eficacia")
report.append("")
report.append(
    "Se realizó una búsqueda dirigida de los datasets "
    "Horlbeck2016, Gasperini2019 y Replogle2022."
)
report.append("")
report.append(
    "**No se realizó ningún merge ni modificación de la matriz Phase 1.**"
)
report.append("")

for dataset in DATASETS:

    report.append(f"## {dataset}")
    report.append("")

    subset = [
        r for r in results
        if r["dataset"] == dataset
    ]

    if not subset:
        report.append("- No se encontraron archivos candidatos.")
        report.append("")
        continue

    for r in subset:

        report.append(f"### `{r['file']}`")
        report.append("")

        report.append(
            f"- Candidatos de score: "
            f"`{r['score_candidates'] or 'NONE'}`"
        )

        report.append(
            f"- Candidatos de llave: "
            f"`{r['key_candidates'] or 'NONE'}`"
        )

        report.append("")

report.append("## Estado")
report.append("")
report.append(
    "- Sanson2018 ya posee `score_norm` en "
    "`data/interim/sanson2018_recovered_official.csv`."
)
report.append(
    "- Horlbeck2016, Gasperini2019 y Replogle2022 "
    "fueron inspeccionados en busca de sus variables originales."
)
report.append(
    "- Todavía no se generó `efficacy_score`."
)
report.append(
    "- Todavía no se generó "
    "`sgRNA_feature_matrix_phase1_with_efficacy.csv`."
)
report.append(
    "- No se realizaron imputaciones, normalizaciones ni merges."
)

with open(REPORT, "w") as fh:
    fh.write("\n".join(report) + "\n")

print("\n" + "=" * 70)
print("REPORT")
print("=" * 70)
print(REPORT)

