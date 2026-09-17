#!/usr/bin/env python3

import pandas as pd


INPUT = "data/interim/features/sgRNA_with_f2_f3_repressive.csv"


print("=" * 60)
print("ChromaCRISPR Phase 1")
print("Week 7 - F5 Genomic Context Diagnosis")
print("=" * 60)


print("\nLoading dataset:")
print(INPUT)

df = pd.read_csv(INPUT)

print("\nDataset loaded")
print(f"Rows    : {len(df):,}")
print(f"Columns : {len(df.columns)}")


print("\nAvailable columns:")
for i, col in enumerate(df.columns):
    print(f"{i:02d}: {col}")


print("\nData types:")
print(df.dtypes)


print("\nFirst rows:")
print(df.head())


# Identify coordinate-like columns

keywords = [
    "chr",
    "chrom",
    "pos",
    "start",
    "end",
    "strand",
    "location",
    "region",
    "coordinate"
]

print("\nPossible genomic context columns:")

for col in df.columns:
    if any(k in col.lower() for k in keywords):
        print("-", col)


print("\nDataset source distribution:")

if "dataset" in df.columns:
    print(df["dataset"].value_counts())

else:
    print("No dataset column detected")


print("\nUnique chromosomes if available:")

for col in df.columns:
    if "chr" in col.lower() or "chrom" in col.lower():
        print(col)
        print(df[col].head())
        print(df[col].nunique())


print("\nDiagnosis completed.")