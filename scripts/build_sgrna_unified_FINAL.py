#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

FILES = [
    ROOT / "data" / "processed" / "sanson2018_reference_FINAL.csv",
    ROOT / "data" / "processed" / "replogle2022_reference_FINAL.csv",
    ROOT / "data" / "processed" / "horlbeck2016_reference_FINAL.csv",
    ROOT / "data" / "processed" / "gasperini2019_reference_FINAL.csv",
]


COMMON_COLUMNS = [
    "dataset",
    "experiment",
    "guide_sequence",
    "gene_symbol",
    "gene_id",
    "chromosome",
    "coordinate",
    "strand",
    "guide_length",
    "score_raw",
    "score_norm",
    "score_type",
    "genome_build",
    "qc_pass",
]


dfs = []

for f in FILES:

    print("\nLoading:", f.name)

    df = pd.read_csv(f)

    # Añadir columnas faltantes
    for col in COMMON_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA

    # Metadatos para Sanson/Replogle
    if df["dataset"].isna().all():

        name = f.stem.replace(
            "_reference_FINAL",
            ""
        )

        df["dataset"] = name.replace(
            "2018",
            "2018"
        ).replace(
            "2022",
            "2022"
        )

    if df["genome_build"].isna().all():
        df["genome_build"] = "hg38"

    if df["qc_pass"].isna().all():
        df["qc_pass"] = True

    print("Rows:", len(df))
    print("Build:", df["genome_build"].unique())

    dfs.append(
        df[COMMON_COLUMNS]
    )


unified = pd.concat(
    dfs,
    ignore_index=True
)


print("\n=== Unified reference ===")
print("Rows:", len(unified))

print("\nDatasets:")
print(
    unified["dataset"]
    .value_counts()
)

print("\nGenome build:")
print(
    unified["genome_build"]
    .value_counts()
)

print("\nDuplicados guía+gen+dataset:")
print(
    unified.duplicated(
        [
            "guide_sequence",
            "gene_symbol",
            "dataset"
        ],
        keep=False
    ).sum()
)


output = (
    ROOT /
    "data" /
    "processed" /
    "sgRNA_unified_FINAL.csv"
)


unified.to_csv(
    output,
    index=False
)


print("\nSaved:", output)