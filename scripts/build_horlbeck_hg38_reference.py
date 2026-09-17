#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(ROOT))

from scripts.parsers.horlbeck import load

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


bed_file = ROOT / "data" / "interim" / "horlbeck.hg38.bed"


# Cargar metadata original hg19
df = load()


# Cargar coordenadas nuevas hg38
bed = pd.read_csv(
    bed_file,
    sep="\t",
    header=None,
    names=[
        "chromosome",
        "start",
        "end",
        "id"
    ]
)


bed["coordinate"] = bed["start"] + 1


# recuperar índice original
bed["index"] = (
    bed["id"]
    .str.split("|")
    .str[-1]
    .astype(int)
)


# ordenar y unir
bed = bed.sort_values("index")

df = df.reset_index(drop=True)

df["chromosome"] = bed["chromosome"].values
df["coordinate"] = bed["coordinate"].values

df["genome_build"] = "hg38"


output = ROOT / "data" / "processed" / "horlbeck2016_reference_FINAL.csv"

df[COMMON_COLUMNS].to_csv(
    output,
    index=False
)

print(output)
print("Rows:", len(df))
print(df["genome_build"].value_counts())