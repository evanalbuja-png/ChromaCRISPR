#!/usr/bin/env python3

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.parsers.gasperini import load

df = load()

print("Filas:", len(df))
print("Coordinate NA:", df["coordinate"].isna().sum())
print("Chromosome NA:", df["chromosome"].isna().sum())

print("\nPrimeras filas con coordinate NA:")
print(
    df[df["coordinate"].isna()][
        ["experiment", "guide_sequence", "chromosome", "coordinate"]
    ].head(20)
)

gasperini = load()

print("Filas:", len(gasperini))
print("Coordinate NA:", gasperini["coordinate"].isna().sum())
print("Chromosome NA:", gasperini["chromosome"].isna().sum())

gasperini = gasperini.dropna(
    subset=["chromosome", "coordinate"]
).copy()

bed = gasperini.copy()

bed["start"] = bed["coordinate"].astype(int) - 1
bed["end"] = bed["coordinate"].astype(int)

bed["name"] = (
    bed["guide_sequence"]
    + "|"
    + bed.index.astype(str).astype(str)
)

bed = bed[
    [
        "chromosome",
        "start",
        "end",
        "name",
    ]
]

output = ROOT / "data/interim/gasperini.hg19.bed"

bed.to_csv(
    output,
    sep="\t",
    header=False,
    index=False,
)

print(output)
print("Rows:", len(bed))