import os
import pandas as pd


F2_FILE = "data/interim/features/sgRNA_with_f2.csv"
F3_FILE = "data/interim/features/f3_h3k27ac_features.csv"
OUT = "data/interim/features/sgRNA_with_f2_f3.csv"


def main():

    print("Loading F2 dataset...")
    guides = pd.read_csv(
        F2_FILE,
        low_memory=False
    )

    print("Loading H3K27ac features...")
    h3k27ac = pd.read_csv(
        F3_FILE,
        low_memory=False
    )


    print("\nCreating region_id for guides...")

    guides["coordinate"] = pd.to_numeric(
        guides["coordinate"],
        errors="coerce"
    )

    valid = guides["coordinate"].notna()

    excluded = (~valid).sum()

    print(
        "Guides without valid coordinate excluded:",
        excluded
    )

    guides = guides[valid].copy()

    guides["coordinate"] = (
        guides["coordinate"]
        .astype(int)
    )


    # Same correction used during F2 integration:
    # CSV coordinate is 1-based, BED/deepTools is 0-based
    guides["start_0based"] = (
        guides["coordinate"] - 1
    )

    guides["end_0based"] = (
        guides["start_0based"]
        + 1
    )

    guides["region_id"] = (
        guides["chromosome"]
        + ":"
        + guides["start_0based"].astype(str)
        + "-"
        + guides["end_0based"].astype(str)
    )


    print("\nGuide regions:")
    print(len(guides))


    print("\nChecking H3K27ac duplicated region_id...")

    duplicates = (
        h3k27ac["region_id"]
        .duplicated()
        .sum()
    )

    print(
        "H3K27ac duplicated regions:",
        duplicates
    )

    h3k27ac = (
        h3k27ac
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    print(
        "Unique H3K27ac regions:",
        len(h3k27ac)
    )


    print("\nJoining H3K27ac features...")

    merged = guides.merge(
        h3k27ac[
            [
                "region_id",
                "H3K27ac_mean",
                "H3K27ac_max",
                "H3K27ac_p90",
                "H3K27ac_sum"
            ]
        ],
        on="region_id",
        how="left",
        validate="many_to_one"
    )


    matched = (
        merged["H3K27ac_mean"]
        .notna()
        .sum()
    )

    coverage = (
        matched /
        len(merged)
        * 100
    )


    print("\nIntegration summary:")
    print("Final guides:", len(merged))
    print(
        "Guides with H3K27ac:",
        matched
    )
    print(
        f"H3K27ac coverage: {coverage:.3f}%"
    )


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
