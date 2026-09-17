#!/usr/bin/env python3
"""
Inspect CRISPR supplementary datasets.

Genera un reporte con:
- hojas disponibles
- dimensiones
- columnas
- tipos de datos
- porcentaje de NA
- primeras filas

Salida:
    logs/dataset_schema_report.md
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data" / "raw" / "crispr_datasets"
LOG_DIR = ROOT / "logs"

LOG_DIR.mkdir(exist_ok=True)

REPORT = LOG_DIR / "dataset_schema_report.md"

xlsx_files = sorted(DATA_DIR.rglob("*.xlsx"))

with open(REPORT, "w", encoding="utf-8") as out:

    out.write("# CRISPR Dataset Schema Report\n\n")

    for file in xlsx_files:

        out.write(f"# {file.name}\n\n")

        try:
            excel = pd.ExcelFile(file)

        except Exception as e:
            out.write(f"ERROR: {e}\n\n")
            continue

        out.write(f"Path: `{file}`\n\n")
        out.write(f"Sheets ({len(excel.sheet_names)}):\n")

        for s in excel.sheet_names:
            out.write(f"- {s}\n")

        out.write("\n---\n")

        for sheet in excel.sheet_names:

            out.write(f"## Sheet: {sheet}\n\n")

            try:
                df = pd.read_excel(file, sheet_name=sheet)

            except Exception as e:
                out.write(f"Cannot read sheet: {e}\n\n")
                continue

            out.write(f"Rows: {len(df):,}\n\n")
            out.write(f"Columns: {len(df.columns)}\n\n")

            out.write("|Column|dtype|Missing %|\n")
            out.write("|---|---|---|\n")

            for c in df.columns:

                miss = df[c].isna().mean() * 100

                out.write(
                    f"|{c}|{df[c].dtype}|{miss:.2f}|\n"
                )

            out.write("\n### First five rows\n\n")

            out.write(df.head().to_markdown(index=False))

            out.write("\n\n")

print(f"\nReport written to:\n{REPORT}")