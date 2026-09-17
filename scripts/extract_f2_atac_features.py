#!/usr/bin/env python3

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


ATAC_MATRIX = Path("data/interim/features/atac_matrix.gz")
BED_FILE = Path("data/interim/features/sgRNA_hg38.bed")

OUTPUT = Path("data/interim/features/f2_atac_features.csv")


def read_compute_matrix(path):

    print(f"Reading deepTools matrix: {path}")

    with gzip.open(path, "rt") as f:

        # deepTools añade '@' antes del JSON
        header = f.readline().strip()

        if header.startswith("@"):
            header = header[1:]

        metadata = json.loads(header)

        print("\nDeepTools parameters:")
        for key, value in metadata.items():
            print(f"{key}: {value}")

        rows = []

        for line in f:
            if line.strip():
                rows.append(
                    line.rstrip("\n").split("\t")
                )

    return metadata, rows


def parse_matrix_rows(rows):

    print("\nParsing matrix rows...")

    parsed = []

    for row in rows:

        # deepTools format:
        # chr start end region score strand signal1...signalN

        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        region = row[3]

        signal = np.array(
            [
                float(x)
                for x in row[6:]
                if x != ""
            ],
            dtype=float
        )

        parsed.append(
            {
                "chromosome": chrom,
                "start": start,
                "end": end,
                "region": region,
                "ATAC_mean": np.mean(signal),
                "ATAC_max": np.max(signal),
                "ATAC_p90": np.percentile(signal, 90),
                "ATAC_sum": np.sum(signal),
                "n_bins": len(signal)
            }
        )

    return pd.DataFrame(parsed)


def load_bed():

    print("\nLoading BED...")

    bed = pd.read_csv(
        BED_FILE,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "start",
            "end",
            "guide_sequence"
        ]
    )

    return bed


def associate_with_bed(atac_df, bed):

    print("\nAssociating ATAC features with BED using region_id...")

    # Crear identificadores únicos
    atac_df["region_id"] = (
        atac_df["chromosome"].astype(str)
        + ":"
        + atac_df["start"].astype(str)
        + "-"
        + atac_df["end"].astype(str)
    )

    bed["region_id"] = (
        bed["chromosome"].astype(str)
        + ":"
        + bed["start"].astype(str)
        + "-"
        + bed["end"].astype(str)
    )

    print("\nChecking duplicated region_id...")

    atac_duplicates = (
        atac_df["region_id"].duplicated().sum()
    )

    bed_duplicates = (
        bed["region_id"].duplicated().sum()
    )

    print(
        f"ATAC duplicated region_id: {atac_duplicates}"
    )

    print(
        f"BED duplicated region_id: {bed_duplicates}"
    )

    # Mantener primera ocurrencia
    atac_df = (
        atac_df
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    bed = (
        bed
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    print("\nAfter deduplication:")
    print(
        f"ATAC unique regions: {len(atac_df)}"
    )
    print(
        f"BED unique regions: {len(bed)}"
    )

    # Inner join 1:1
    merged = bed.merge(
        atac_df,
        on="region_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_bed", "_atac")
    )

    print("\nFinal association:")
    print(
        f"Matched regions: {len(merged)}"
    )

    discarded_bed = len(bed) - len(merged)

    print(
        f"BED regions discarded: {discarded_bed}"
    )

    print(
        "Coverage:",
        f"{len(merged)/len(bed)*100:.3f}%"
    )

    return merged


def main():

    metadata, rows = read_compute_matrix(
        ATAC_MATRIX
    )

    print(
        f"\nMatrix rows loaded: {len(rows)}"
    )

    atac_features = parse_matrix_rows(rows)

    print(
        f"Regions processed: {len(atac_features)}"
    )

    bed = load_bed()

    print(
        f"BED regions: {len(bed)}"
    )

    final = associate_with_bed(
        atac_features,
        bed
    )

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    final.to_csv(
        OUTPUT,
        index=False
    )

    print("\nSaved:")
    print(OUTPUT)

    print("\nCoverage:")
    print(
        f"Regions with ATAC features: "
        f"{len(final)}/{len(bed)}"
    )


if __name__ == "__main__":
    main()
