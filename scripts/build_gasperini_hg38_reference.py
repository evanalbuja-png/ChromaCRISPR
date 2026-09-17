#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(ROOT))

from scripts.parsers.gasperini import load


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


# Metadata original hg19
df = load()


# BED generado para liftOver
bed_file = (
    ROOT /
    "data" /
    "interim" /
    "gasperini.hg38.bed"
)


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


# Recuperar índice original dentro del parser
bed["index"] = (
    bed["id"]
    .str.split("|")
    .str[-1]
    .astype(int)
)


bed = bed.sort_values("index")


# Mantener solo guías con coordenada válida
df = df.reset_index(drop=True)

df = df.loc[
    bed["index"]
].copy()


df["chromosome"] = bed["chromosome"].values
df["coordinate"] = bed["coordinate"].values

df["genome_build"] = "hg38"


output = (
    ROOT /
    "data" /
    "processed" /
    "gasperini2019_reference_FINAL.csv"
)


df[COMMON_COLUMNS].to_csv(
    output,
    index=False
)


print(output)
print("Rows:", len(df))
print(df["genome_build"].value_counts())