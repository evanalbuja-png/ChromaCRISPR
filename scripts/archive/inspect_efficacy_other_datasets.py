import pandas as pd
from pathlib import Path
import re
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

BASE = Path("data/raw/crispr_datasets")

DATASETS = {
    "Sanson2018": BASE / "Sanson2018",
    "Replogle2022": BASE / "Replogle2022",
    "Gasperini2019": BASE / "Gasperini2019",
}

EFFICACY_PATTERNS = re.compile(
    r"(effic|activ|activity|phenotype|phenotyp|"
    r"log2fc|logfc|lfc|depletion|enrichment|"
    r"fold.?change|fitness|effect|effect_size|"
    r"score|beta|z.?score|growth|abundance)",
    re.I,
)

JOIN_PATTERNS = re.compile(
    r"(guide|sgrna|grna|sequence|seq|"
    r"gene|symbol|target|ensembl|"
    r"sgid|spacer)",
    re.I,
)


def detect_header_rows(path, sheet, nrows=20):
    """
    Busca filas que parezcan contener nombres de columnas reales.
    """
    raw = pd.read_excel(
        path,
        sheet_name=sheet,
        header=None,
        nrows=nrows,
    )

    candidates = []

    for i, row in raw.iterrows():
        values = [str(x).strip() for x in row.tolist() if pd.notna(x)]

        score = 0

        for value in values:
            if re.search(
                r"(guide|sgrna|grna|sequence|gene|symbol|"
                r"count|read|pDNA|HT29|A375|log|fold|score|"
                r"depletion|enrichment)",
                value,
                re.I,
            ):
                score += 1

        if score >= 2:
            candidates.append((i, score, values))

    return candidates


def read_excel_sheet(path, sheet):
    """
    Lee una hoja intentando identificar automáticamente el header real.
    """
    candidates = detect_header_rows(path, sheet)

    if candidates:
        header_row = max(candidates, key=lambda x: x[1])[0]
    else:
        header_row = 0

    df = pd.read_excel(
        path,
        sheet_name=sheet,
        header=header_row,
    )

    # Eliminar columnas completamente vacías
    df = df.dropna(axis=1, how="all")

    # Normalizar nombres
    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    return df, header_row, candidates


def inspect_dataframe(path, sheet, df, header_row):
    print("\n" + "=" * 100)
    print(f"FILE: {path}")
    print(f"SHEET: {sheet}")
    print(f"HEADER ROW: {header_row}")
    print(f"ROWS: {len(df)}")
    print("=" * 100)

    print("\nCOLUMNS:")
    for c in df.columns:
        print(f"  {c!r}")

    efficacy = [
        c for c in df.columns
        if EFFICACY_PATTERNS.search(str(c))
    ]

    numeric = [
        c for c in df.columns
        if pd.api.types.is_numeric_dtype(df[c])
    ]

    joins = [
        c for c in df.columns
        if JOIN_PATTERNS.search(str(c))
    ]

    print("\n--- EFFICACY / ACTIVITY CANDIDATES ---")

    if not efficacy:
        print("  NONE DETECTED")
    else:
        for c in efficacy:
            s = df[c]

            print(f"\n  COLUMN: {c!r}")
            print(f"    dtype: {s.dtype}")
            print(f"    numeric: {pd.api.types.is_numeric_dtype(s)}")
            print(f"    non_null: {s.notna().sum()}")
            print(f"    unique: {s.nunique(dropna=True)}")

            if pd.api.types.is_numeric_dtype(s):
                print(f"    min: {s.min()}")
                print(f"    max: {s.max()}")
                print(f"    mean: {s.mean()}")
                print(f"    median: {s.median()}")

    print("\n--- ALL NUMERIC COLUMNS ---")

    if not numeric:
        print("  NONE")
    else:
        for c in numeric:
            s = df[c]

            print(
                f"  {c!r} | "
                f"non_null={s.notna().sum()} | "
                f"unique={s.nunique(dropna=True)} | "
                f"min={s.min()} | "
                f"max={s.max()} | "
                f"mean={s.mean()}"
            )

    print("\n--- POSSIBLE JOIN COLUMNS ---")

    if not joins:
        print("  NONE DETECTED")
    else:
        for c in joins:
            s = df[c]

            print(
                f"  {c!r} | "
                f"dtype={s.dtype} | "
                f"non_null={s.notna().sum()} | "
                f"unique={s.nunique(dropna=True)}"
            )

    print("\n--- FIRST 5 ROWS ---")
    print(df.head(5).to_string(index=False))


def inspect_dataset(name, directory):

    print("\n\n")
    print("#" * 100)
    print(f"# DATASET: {name}")
    print(f"# DIRECTORY: {directory}")
    print("#" * 100)

    files = sorted(
        f for f in directory.rglob("*")
        if f.is_file()
        and f.suffix.lower() in {".xlsx", ".xls"}
    )

    for path in files:

        print("\n")
        print(f"### FILE: {path}")

        try:
            xls = pd.ExcelFile(path)

            print("SHEETS:")
            for sheet in xls.sheet_names:
                print(f"  - {sheet}")

            for sheet in xls.sheet_names:

                try:
                    df, header_row, candidates = read_excel_sheet(
                        path,
                        sheet,
                    )

                    inspect_dataframe(
                        path,
                        sheet,
                        df,
                        header_row,
                    )

                except Exception as e:
                    print(
                        f"\nERROR reading sheet "
                        f"{sheet!r}: {e}"
                    )

        except Exception as e:
            print(f"ERROR opening {path}: {e}")


def main():

    print("=" * 100)
    print("CHROMACRISPR PHASE 1")
    print("STEP 2 — EFFICACY INSPECTION")
    print("=" * 100)

    for name, directory in DATASETS.items():

        if not directory.exists():
            print(
                f"\nWARNING: dataset directory does not exist: "
                f"{directory}"
            )
            continue

        inspect_dataset(
            name,
            directory,
        )

    print("\n")
    print("=" * 100)
    print("STEP 2 INSPECTION FINISHED")
    print("=" * 100)


if __name__ == "__main__":
    main()
