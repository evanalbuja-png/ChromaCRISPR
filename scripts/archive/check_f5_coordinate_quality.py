#!/usr/bin/env python3

import pandas as pd


INPUT = "data/interim/features/sgRNA_with_f2_f3_repressive.csv"


df = pd.read_csv(
    INPUT,
    low_memory=False
)


print("="*60)
print("F5 Coordinate Quality Check")
print("="*60)


print("\nGenome build distribution:")
print(
    df["genome_build"]
    .value_counts(dropna=False)
)


print("\nMissing coordinate information:")

cols = [
    "chromosome",
    "coordinate",
    "start",
    "end",
    "strand"
]

for col in cols:
    missing = df[col].isna().sum()
    print(
        f"{col}: {missing}"
    )


print("\nCoordinate ranges:")

print(df["coordinate"].describe())

print("\nStart:")
print(df["start"].describe())

print("\nEnd:")
print(df["end"].describe())


print("\nChromosome distribution:")
print(
    df["chromosome"]
    .value_counts()
    .head(30)
)


print("\nStrand distribution:")
print(
    df["strand"]
    .value_counts(dropna=False)
)


print("\nDataset vs genome build:")
print(
    pd.crosstab(
        df["dataset"],
        df["genome_build"]
    )
)


print("\nCompleted.")
