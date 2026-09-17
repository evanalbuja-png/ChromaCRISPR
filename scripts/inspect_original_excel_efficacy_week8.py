#!/usr/bin/env python3

import os
import pandas as pd

FILES = {
    "Horlbeck2016": "data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx",
    "Sanson2018": "data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx",
    "Gasperini2019_S1": "data/raw/crispr_datasets/Gasperini2019/NIHMS1038673-supplement-TableS1.xlsx",
    "Gasperini2019_S2": "data/raw/crispr_datasets/Gasperini2019/NIHMS1038673-supplement-TableS2.xlsx",
    "Replogle2022": "data/raw/crispr_datasets/Replogle2022/NIHMS1812939-supplement-11.xlsx",
}

REPORT = "logs/week8_original_excel_inspection.md"

score_terms = [
    "score", "efficacy", "efficiency", "activity",
    "log2fc", "log2_fc", "phenotype", "effect",
    "lfc", "fold", "depletion", "enrichment"
]

guide_terms = [
    "guide", "sgrna", "sgRNA", "sequence", "spacer",
    "protospacer", "grna"
]

def matches(columns, terms):
    result = []
    for col in columns:
        s = str(col).lower()
        if any(term.lower() in s for term in terms):
            result.append(str(col))
    return result


lines = []

lines.append("# ChromaCRISPR Phase 1 — Week 8")
lines.append("")
lines.append("## Inspección de datasets originales de eficacia")
lines.append("")

for name, path in FILES.items():

    lines.append(f"## {name}")
    lines.append("")
    lines.append(f"Archivo: `{path}`")
    lines.append("")

    print("=" * 70)
    print(name)
    print(path)
    print("=" * 70)

    if not os.path.exists(path):
        print("ERROR: archivo no encontrado")
        lines.append("**ERROR: archivo no encontrado.**")
        lines.append("")
        continue

    try:
        xls = pd.ExcelFile(path)
    except Exception as e:
        print("ERROR:", e)
        lines.append(f"**ERROR al abrir Excel:** `{e}`")
        lines.append("")
        continue

    print("Hojas:")
    for sheet in xls.sheet_names:
        print(" -", sheet)

    lines.append("### Hojas")
    lines.append("")
    for sheet in xls.sheet_names:
        lines.append(f"- `{sheet}`")
    lines.append("")

    for sheet in xls.sheet_names:

        print()
        print("--- SHEET:", sheet, "---")

        try:
            df = pd.read_excel(
                path,
                sheet_name=sheet,
                nrows=8
            )
        except Exception as e:
            print("READ ERROR:", e)
            lines.append(f"### Hoja `{sheet}`")
            lines.append("")
            lines.append(f"READ ERROR: `{e}`")
            lines.append("")
            continue

        columns = list(df.columns)

        guide_cols = matches(columns, guide_terms)
        score_cols = matches(columns, score_terms)

        print("Columns:")
        for i, col in enumerate(columns):
            print(f"  {i}: {col}")

        print("Potential GUIDE columns:")
        for col in guide_cols:
            print("  -", col)

        print("Potential SCORE/EFFICACY columns:")
        for col in score_cols:
            print("  -", col)

        print("\nPreview:")
        print(df.head(3).to_string(index=False))

        lines.append(f"### Hoja `{sheet}`")
        lines.append("")
        lines.append("**Columnas:**")
        lines.append("")

        for i, col in enumerate(columns):
            lines.append(f"- `{i}`: `{col}`")

        lines.append("")
        lines.append("**Posibles columnas de guía:**")
        lines.append("")

        if guide_cols:
            for col in guide_cols:
                lines.append(f"- `{col}`")
        else:
            lines.append("- Ninguna detectada automáticamente.")

        lines.append("")
        lines.append("**Posibles columnas de eficacia/score:**")
        lines.append("")

        if score_cols:
            for col in score_cols:
                lines.append(f"- `{col}`")
        else:
            lines.append("- Ninguna detectada automáticamente.")

        lines.append("")
        lines.append("**Primeras filas:**")
        lines.append("")
        lines.append("```text")
        lines.append(df.head(3).to_string(index=False))
        lines.append("```")
        lines.append("")

print()
print("=" * 70)
print("REPORT")
print("=" * 70)

with open(REPORT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(REPORT)
