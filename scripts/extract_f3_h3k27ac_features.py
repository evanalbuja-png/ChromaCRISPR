import gzip
import json
import os
import numpy as np
import pandas as pd


MATRIX = "data/interim/features/histones/H3K27ac_matrix.gz"
BED = "data/interim/features/sgRNA_hg38.bed"
OUT = "data/interim/features/f3_h3k27ac_features.csv"


def read_matrix(path):
    print("Reading deepTools matrix:", path)

    with gzip.open(path, "rt") as f:
        header = f.readline().strip()

        if header.startswith("@"):
            header = header[1:]

        metadata = json.loads(header)

        rows = []
        for line in f:
            if line.strip():
                rows.append(line.strip().split())

    return metadata, rows


def load_bed(path):
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chromosome", "start", "end", "guide_sequence"]
    )

    bed["region_id"] = (
        bed["chromosome"]
        + ":"
        + bed["start"].astype(str)
        + "-"
        + bed["end"].astype(str)
    )

    return bed


def extract_features(rows):
    records = []

    for row in rows:
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])

        region_id = f"{chrom}:{start}-{end}"

        signal = np.array(
            row[6:],
            dtype=float
        )

        records.append(
            {
                "region_id": region_id,
                "chromosome": chrom,
                "start": start,
                "end": end,
                "H3K27ac_mean": np.mean(signal),
                "H3K27ac_max": np.max(signal),
                "H3K27ac_p90": np.percentile(signal, 90),
                "H3K27ac_sum": np.sum(signal),
                "n_bins": len(signal)
            }
        )

    return pd.DataFrame(records)


def main():

    metadata, rows = read_matrix(MATRIX)

    print("\nDeepTools parameters:")
    for k, v in metadata.items():
        print(f"{k}: {v}")

    print("\nMatrix rows loaded:", len(rows))

    print("\nExtracting H3K27ac features...")
    features = extract_features(rows)

    print("Regions extracted:", len(features))

    print("\nChecking duplicated region_id...")

    duplicated = features["region_id"].duplicated().sum()

    print("Duplicated region_id:", duplicated)

    features = (
        features
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    print("Unique matrix regions:", len(features))


    print("\nLoading BED...")
    bed = load_bed(BED)

    print("BED regions:", len(bed))

    bed_unique = (
        bed
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    print("Unique BED regions:", len(bed_unique))


    print("\nAssociating H3K27ac with BED...")

    merged = bed_unique.merge(
        features,
        on="region_id",
        how="inner",
        suffixes=("_bed", "_matrix")
    )


    coverage = (
        len(merged) /
        len(bed_unique)
        * 100
    )

    print("\nFinal association:")
    print("Matched regions:", len(merged))
    print("BED discarded:", len(bed_unique)-len(merged))
    print(f"Coverage: {coverage:.3f}%")


    os.makedirs(
        os.path.dirname(OUT),
        exist_ok=True
    )

    merged.to_csv(
        OUT,
        index=False
    )

    print("\nSaved:")
    print(OUT)


if __name__ == "__main__":
    main()
