#!/usr/bin/env python3

import os
import pandas as pd

BASE = "data/raw/crispr_datasets"
REPORT = "logs/week8_efficacy_extraction_inspection.md"

DATASETS = {
    "Horlbeck2016": os.path.join(
        BASE, "Horlbeck2016", "elife-19760-supp1-v2.xlsx"
    ),
    "Sanson2018": os.path.join(
        BASE, "Sanson2018", "41467_2018_7901_MOESM6_ESM.xlsx"
    ),
    "Gasperini2019_S1": os.path.join(
        BASE, "Gasperini2019", "NIHMS1038673-supplement-TableS1.xlsx"
    ),
    "Gasperini2019_S2": os.path.join(
        BASE, "Gasperini2019", "NIHMS1038673-supplement-TableS2.xlsx"
    ),
    "Replogle2022": os.path.join(
        BASE, "Replogle2022", "NIHMS1812939-supplement-11.xlsx"
    ),
}

def inspect_excel(name, path):
    print("\n" + "=" * 80)
    print(name)
    print("=" * 80)

    if not os.path.exists(path):
        print("FILE NOT FOUND:", path)
        return []

    xls = pd.ExcelFile(path)
    print("Sheets:")
    for s in xls.sheet_names:
        print("  -", s)

    results = []

    for sheet in xls.sheet_names:
        print("\n---", sheet, "---")

        try:
            df = pd.read_excel(path, sheet_name=sheet, header=None)
        except Exception as e:
            print("READ ERROR:", e)
            continue

        print("Shape:", df.shape)

        # Mostrar primeras 8 filas para entender la estructura real
        print("\nFirst rows:")
        print(df.head(8).to_string(index=False, header=False))

        # Buscar texto relevante en todo el contenido
        keywords = [
            "score",
            "activity",
            "efficacy",
            "phenotype",
            "fold",
            "log2",
            "log",
            "depletion",
            "enrichment",
            "count",
            "read",
            "fitness",
            "effect",
            "sgRNA",
            "Spacer",
            "targeting sequence",
            "Target_Site",
        ]

        matches = []

        for i in range(min(len(df), 100)):
            for j in range(df.shape[1]):
                value = df.iat[i, j]

                if pd.isna(value):
                    continue

                text = str(value).strip()

                if any(k.lower() in text.lower() for k in keywords):
                    matches.append((i, j, text))

        if matches:
            print("\nRelevant cells:")
            for i, j, text in matches[:100]:
                print(f"  row={i}, col={j}: {text}")

        results.append({
            "dataset": name,
            "sheet": sheet,
            "rows": df.shape[0],
            "columns": df.shape[1],
            "matches": matches,
        })

    return results


all_results = []

for name, path in DATASETS.items():
    all_results.extend(inspect_excel(name, path))


# ----------------------------------------------------------------------
# Markdown report
# ----------------------------------------------------------------------

os.makedirs(os.path.dirname(REPORT), exist_ok=True)

with open(REPORT, "w", encoding="utf-8") as f:

    f.write("# ChromaCRISPR Phase 1 — Week 8\n\n")
    f.write("## Inspección para extracción de eficacia desde datasets originales\n\n")

    for r in all_results:

        f.write(f"## {r['dataset']} — `{r['sheet']}`\n\n")
        f.write(f"- Filas: {r['rows']}\n")
        f.write(f"- Columnas: {r['columns']}\n\n")

        if r["matches"]:
            f.write("### Celdas relevantes detectadas\n\n")

            for i, j, text in r["matches"][:100]:
                safe = text.replace("|", "\\|")
                f.write(f"- row `{i}`, col `{j}`: `{safe}`\n")

        else:
            f.write("No se detectaron términos relevantes automáticamente.\n")

        f.write("\n")

print("\n" + "=" * 80)
print("REPORT SAVED:", REPORT)
print("=" * 80)

