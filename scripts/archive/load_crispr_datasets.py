#!/usr/bin/env python3
"""
ChromaCRISPR Phase 1

Load and QC CRISPR datasets
"""

from pathlib import Path

import pandas as pd

from parsers.horlbeck import load as load_horlbeck
from parsers.sanson import load as load_sanson
from parsers.gasperini import load as load_gasperini
from parsers.replogle import load as load_replogle

ROOT = Path(__file__).resolve().parents[1]

OUTDIR = ROOT / "data" / "interim"
OUTDIR.mkdir(exist_ok=True)


def qc(df):
    """
    Control de calidad de las sgRNAs.

    Aplica los filtros en el siguiente orden:
    1. Secuencias faltantes.
    2. Longitud mínima (19 nt).
    3. Caracteres válidos (A/C/G/T).
    4. Duplicados.

    Imprime un reporte resumido y devuelve el DataFrame filtrado.
    """

    initial = len(df)

    print("\n==========================")
    print("QUALITY CONTROL REPORT")
    print("==========================")

    # ------------------------------------------------------------------
    # 1. Secuencias faltantes
    # ------------------------------------------------------------------

    missing_mask = df["guide_sequence"].isna()
    missing = missing_mask.sum()

    df = df[~missing_mask]

    # ------------------------------------------------------------------
    # 2. Longitud mínima
    # ------------------------------------------------------------------

    length_mask = df["guide_sequence"].str.len() >= 19
    invalid_length = (~length_mask).sum()

    df = df[length_mask]

    # ------------------------------------------------------------------
    # 3. Solo A,C,G,T
    # ------------------------------------------------------------------

    sequence_mask = df["guide_sequence"].str.fullmatch(
        r"[ACGT]+",
        na=False,
    )

    invalid_sequence = (~sequence_mask).sum()

    df = df[sequence_mask]

    # ------------------------------------------------------------------
    # 4. Duplicados
    # ------------------------------------------------------------------

    duplicates = df.duplicated(
        subset="guide_sequence"
    ).sum()

    df.drop_duplicates(
    subset=[
        "dataset",
        "experiment",
        "guide_sequence",
        "gene_id",
        "gene_symbol"
    ]
    )

    # ------------------------------------------------------------------
    # Resumen
    # ------------------------------------------------------------------

    final = len(df)
    removed = initial - final
    retention = (final / initial) * 100

    print(f"Initial guides:        {initial:,}")
    print(f"Missing sequence:      {missing:,}")
    print(f"Invalid length:        {invalid_length:,}")
    print(f"Invalid characters:    {invalid_sequence:,}")
    print(f"Duplicate guides:      {duplicates:,}")
    print("--------------------------------")
    print(f"Removed total:         {removed:,}")
    print(f"Final guides:          {final:,}")
    print(f"Retention:             {retention:.2f}%")
    print()

    qc_stats = {
    "initial": int(initial),
    "missing": int(missing),
    "invalid_length": int(invalid_length),
    "invalid_characters": int(invalid_sequence),
    "duplicates": int(duplicates),
    "removed": int(removed),
    "final": int(final),
    "retention": round(float(retention), 2),
    }

    return df, qc_stats

def main():

    print("Loading Horlbeck...")
    horlbeck = load_horlbeck()

    print("Loading Sanson...")
    sanson = load_sanson()

    print("Loading Gasperini...")
    gasperini = load_gasperini()

    print("Loading Replogle...")
    replogle = load_replogle()

    df = pd.concat(
    [
        horlbeck,
        sanson,
        gasperini,
        replogle,
    ],
    ignore_index=True,
    )

    df, qc_stats = qc(df)
    print("\nQC statistics dictionary:")
    print(qc_stats)

    outfile = OUTDIR / "sgRNA_unified.csv"

    dup = (
        df.groupby(
            ["dataset","experiment","guide_sequence"]
        )
        .size()
    )

    print(
        dup[dup>1].head(20)
    )

    df.to_csv(outfile, index=False)

    print()

    print(df.head())

    print()

    print(df.describe(include="all"))

    print()

    print(f"Saved to {outfile}")


if __name__ == "__main__":
    main()