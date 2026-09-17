import pandas as pd
import numpy as np
from pathlib import Path

FILE = Path(
    "data/raw/crispr_datasets/Sanson2018/"
    "41467_2018_7901_MOESM6_ESM.xlsx"
)

OUT = Path("results/efficacy_inspection")
OUT.mkdir(parents=True, exist_ok=True)

PSEUDOCOUNT = 1


def load_raw(sheet):
    raw = pd.read_excel(FILE, sheet_name=sheet, header=None)

    # Detectar la fila que contiene sgRNA Sequence
    seq_row = raw.index[
        raw.iloc[:, 0].astype(str).str.strip().eq("sgRNA Sequence")
    ][0]

    # Fila inmediatamente anterior:
    # pDNA / HT29 / HT29 / HT29 / A375 / A375 / A375
    condition_row = seq_row - 1

    # Construir nombres explícitos
    columns = ["guide_sequence"]

    for i in range(1, raw.shape[1]):
        condition = str(raw.iloc[condition_row, i]).strip()
        replicate = str(raw.iloc[seq_row, i]).strip()

        if condition == "pDNA":
            name = "pDNA"
        else:
            name = f"{condition}_{replicate}"

        columns.append(name)

    # Datos después de las filas de encabezado
    df = raw.iloc[seq_row + 1:].copy()
    df.columns = columns
    df = df.reset_index(drop=True)

    # Eliminar filas sin secuencia
    df = df[df["guide_sequence"].notna()].copy()

    # Convertir conteos
    for c in df.columns:
        if c != "guide_sequence":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def evaluate(sheet):
    print("\n" + "=" * 100)
    print(f"SHEET: {sheet}")
    print("=" * 100)

    df = load_raw(sheet)

    print(f"ROWS: {len(df)}")
    print(f"COLUMNS: {list(df.columns)}")

    ht29 = [
        "HT29_RepA",
        "HT29_RepB",
        "HT29_RepC"
    ]

    a375 = [
        "A375_RepA",
        "A375_RepB",
        "A375_RepC"
    ]

    print(f"\nHT29 columns: {ht29}")
    print(f"A375 columns: {a375}")

    missing = [
        c for c in ["pDNA"] + ht29 + a375
        if c not in df.columns
    ]

    if missing:
        raise RuntimeError(
            f"Missing expected columns: {missing}"
        )

    # Derivar log2 depletion respecto al input pDNA
    for c in ht29 + a375:
        df[f"log2FC_{c}"] = np.log2(
            (df[c] + PSEUDOCOUNT) /
            (df["pDNA"] + PSEUDOCOUNT)
        )

    ht29_scores = [f"log2FC_{c}" for c in ht29]
    a375_scores = [f"log2FC_{c}" for c in a375]

    df["HT29_log2FC_mean"] = df[ht29_scores].mean(axis=1)
    df["A375_log2FC_mean"] = df[a375_scores].mean(axis=1)

    df["Sanson_log2FC_mean"] = df[
        ["HT29_log2FC_mean", "A375_log2FC_mean"]
    ].mean(axis=1)

    print("\n--- DERIVED SCORE ---")

    for c in [
        "HT29_log2FC_mean",
        "A375_log2FC_mean",
        "Sanson_log2FC_mean"
    ]:
        s = df[c]

        print(
            f"{c}: "
            f"non_null={s.notna().sum()}, "
            f"mean={s.mean():.4f}, "
            f"median={s.median():.4f}, "
            f"min={s.min():.4f}, "
            f"max={s.max():.4f}"
        )

    print("\n--- REPLICATE SPEARMAN CORRELATIONS ---")

    print("\nHT29:")
    print(
        df[ht29_scores]
        .corr(method="spearman")
        .to_string()
    )

    print("\nA375:")
    print(
        df[a375_scores]
        .corr(method="spearman")
        .to_string()
    )

    outfile = (
        OUT /
        f"{sheet.replace(' ', '_')}_derived_scores.tsv"
    )

    df.to_csv(outfile, sep="\t", index=False)

    print(f"\nSAVED: {outfile}")

    return df


setA = evaluate("SetA raw reads")
setB = evaluate("SetB raw reads")

print("\n" + "=" * 100)
print("SANSON2018 — PARSER CORREGIDO")
print("=" * 100)
print("""
Ambos sets contienen:
  pDNA
  HT29 RepA/RepB/RepC
  A375 RepA/RepB/RepC

Se derivó para cada guía:

    log2FC = log2((selected + 1) / (pDNA + 1))

y posteriormente:

    HT29_log2FC_mean
    A375_log2FC_mean
    Sanson_log2FC_mean

El score es un depletion/activity proxy derivado de los raw reads.
""")
