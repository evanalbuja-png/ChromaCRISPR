#!/usr/bin/env python3

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.parsers.horlbeck import load

df = load()

bed = df.copy()

bed["start"] = bed["coordinate"].astype(int) - 1
bed["end"] = bed["coordinate"].astype(int)

bed["name"] = (
    bed["guide_sequence"]
    + "|"
    + bed["gene_symbol"].fillna("NA").astype(str)
    + "|"
    + bed.index.astype(str)
)

bed = bed[
    [
        "chromosome",
        "start",
        "end",
        "name",
    ]
]

output = ROOT / "data/interim/horlbeck.hg19.bed"

bed.to_csv(
    output,
    sep="\t",
    header=False,
    index=False,
)

print(output)
print("Rows:", len(bed))