import gzip
import json
import numpy as np
import pandas as pd


MATRIX = "data/interim/features/histones/H3K27me3_matrix.gz"
BED = "data/interim/features/sgRNA_hg38.bed"
OUT = "data/interim/features/f3_h3k27me3_features.csv"
LOG = "logs/week6_log.md"


def read_matrix(path):

    with gzip.open(path, "rt") as f:

        header = f.readline().strip()

        if header.startswith("@"):
            header = header[1:]

        metadata = json.loads(header)

        rows = []

        for line in f:
            values = line.strip().split()

            if values:
                rows.append(values)

    return metadata, rows


def main():

    print("Reading deepTools matrix...")

    metadata, rows = read_matrix(MATRIX)

    print("\nDeepTools parameters:")

    for k, v in metadata.items():
        print(f"{k}: {v}")

    print(f"\nMatrix rows loaded: {len(rows)}")


    print("\nExtracting H3K27me3 features...")

    features = []

    for row in rows:

        signal = np.array(
            row[-100:],
            dtype=float
        )

        features.append(
            {
                "H3K27me3_mean": np.mean(signal),
                "H3K27me3_max": np.max(signal),
                "H3K27me3_p90": np.percentile(signal, 90),
                "H3K27me3_sum": np.sum(signal),
            }
        )


    atac = pd.DataFrame(features)

    print(
        f"Regions extracted: {len(atac)}"
    )


    print("\nLoading BED...")

    bed = pd.read_csv(
        BED,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "start",
            "end",
            "guide_sequence"
        ],
        dtype={
            "chromosome": str,
            "start": int,
            "end": int
        }
    )


    bed["region_id"] = (
        bed["chromosome"]
        + ":"
        + bed["start"].astype(str)
        + "-"
        + bed["end"].astype(str)
    )


    atac["region_id"] = (
        bed["region_id"]
        .iloc[:len(atac)]
        .values
    )


    print("\nChecking duplicated region_id...")

    duplicated = atac["region_id"].duplicated().sum()

    print(
        f"Duplicated regions: {duplicated}"
    )


    atac = atac.drop_duplicates(
        subset="region_id",
        keep="first"
    )

    print(
        f"Unique H3K27me3 regions: {len(atac)}"
    )


    bed = bed.drop_duplicates(
        subset="region_id",
        keep="first"
    )


    merged = bed.merge(
        atac,
        on="region_id",
        how="inner",
        validate="one_to_one"
    )


    print("\nFinal association:")

    print(
        f"Matched regions: {len(merged)}"
    )

    print(
        f"Coverage: {len(merged)/len(bed)*100:.3f}%"
    )


    merged.to_csv(
        OUT,
        index=False
    )


    print("\nSaved:")
    print(OUT)


    with open(LOG, "a") as f:

        f.write("\n\n## F3 H3K27me3 extraction\n")
        f.write(
            f"- Matrix rows: {len(rows)}\n"
        )
        f.write(
            f"- Unique regions: {len(merged)}\n"
        )
        f.write(
            f"- Coverage: {len(merged)/len(bed)*100:.3f}%\n"
        )
        f.write(
            f"- Duplicate region_id removed: {duplicated}\n"
        )


if __name__ == "__main__":
    main()
