import pandas as pd


F2F3_FILE = "data/interim/features/sgRNA_with_f2_f3.csv"
H3K4_FILE = "data/interim/features/f3_h3k4me3_features.csv"

OUT = "data/interim/features/sgRNA_with_f2_f3_full.csv"
LOG = "logs/week6_log.md"


def main():

    print("Loading F2 + H3K27ac dataset...")
    guides = pd.read_csv(
        F2F3_FILE,
        low_memory=False
    )

    print("Loading H3K4me3 features...")
    h3k4 = pd.read_csv(
        H3K4_FILE,
        low_memory=False
    )


    print("\nCreating region_id for guides...")

    guides["coordinate"] = pd.to_numeric(
        guides["coordinate"],
        errors="coerce"
    )

    invalid = guides["coordinate"].isna().sum()

    print(
        f"Guides without valid coordinate: {invalid}"
    )

    guides = guides[
        guides["coordinate"].notna()
    ].copy()

    guides["coordinate"] = (
        guides["coordinate"]
        .astype(int)
    )


    # CSV coordinates are 1-based
    # BED/deepTools coordinates are 0-based
    guides["start"] = guides["coordinate"] - 1
    guides["end"] = guides["start"] + 1

    guides["region_id"] = (
        guides["chromosome"]
        + ":"
        + guides["start"].astype(str)
        + "-"
        + guides["end"].astype(str)
    )


    print("\nGuide regions:")
    print(len(guides))


    print("\nChecking H3K4me3 duplicated region_id...")

    duplicated = h3k4["region_id"].duplicated().sum()

    print(
        f"H3K4me3 duplicated regions: {duplicated}"
    )

    h3k4 = h3k4.drop_duplicates(
        subset="region_id",
        keep="first"
    )

    print(
        f"Unique H3K4me3 regions: {len(h3k4)}"
    )


    print("\nJoining H3K4me3 features...")


    feature_cols = [
        "region_id",
        "H3K4me3_mean",
        "H3K4me3_max",
        "H3K4me3_p90",
        "H3K4me3_sum"
    ]


    merged = guides.merge(
        h3k4[feature_cols],
        on="region_id",
        how="left",
        validate="many_to_one"
    )


    matched = (
        merged["H3K4me3_mean"]
        .notna()
        .sum()
    )


    print("\nIntegration summary:")
    print(
        f"Final guides: {len(merged)}"
    )

    print(
        f"Guides with H3K4me3: {matched}"
    )

    print(
        f"H3K4me3 coverage: {matched/len(merged)*100:.3f}%"
    )


    merged.to_csv(
        OUT,
        index=False
    )


    print("\nSaved:")
    print(OUT)


    with open(LOG, "a") as f:

        f.write("\n\n## H3K4me3 F3 integration\n")
        f.write(
            f"- Input dataset: {F2F3_FILE}\n"
        )
        f.write(
            f"- H3K4me3 regions: {len(h3k4)}\n"
        )
        f.write(
            f"- Final guides: {len(merged)}\n"
        )
        f.write(
            f"- Guides with H3K4me3: {matched}\n"
        )
        f.write(
            f"- Coverage: {matched/len(merged)*100:.3f}%\n"
        )


if __name__ == "__main__":
    main()
