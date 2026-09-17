import pandas as pd


INPUT = "data/interim/features/sgRNA_with_f2_f3_full.csv"
H3K27ME3 = "data/interim/features/f3_h3k27me3_features.csv"

OUT = "data/interim/features/sgRNA_with_f2_f3_repressive.csv"
LOG = "logs/week6_log.md"


def create_region_id(df):
    df = df.copy()

    df["coordinate"] = pd.to_numeric(
        df["coordinate"],
        errors="coerce"
    )

    valid = df["coordinate"].notna()

    print(
        f"Guides without valid coordinate: {(~valid).sum()}"
    )

    df = df[valid].copy()

    # Corrección 1-based CSV -> 0-based BED/deepTools
    df["start"] = (
        df["coordinate"]
        .astype(int)
        - 1
    )

    df["end"] = df["start"] + 1

    df["region_id"] = (
        df["chromosome"]
        + ":"
        + df["start"].astype(str)
        + "-"
        + df["end"].astype(str)
    )

    return df


def main():

    print("Loading F2 + H3K27ac + H3K4me3 dataset...")
    guides = pd.read_csv(
        INPUT,
        low_memory=False
    )

    print("Loading H3K27me3 features...")
    h3k27me3 = pd.read_csv(
        H3K27ME3,
        low_memory=False
    )

    print("\nCreating region_id for guides...")
    guides = create_region_id(guides)

    print("\nGuide regions:")
    print(len(guides))

    print("\nChecking duplicated H3K27me3 region_id...")

    dup = h3k27me3["region_id"].duplicated().sum()

    print(
        f"H3K27me3 duplicated regions: {dup}"
    )

    if dup > 0:
        h3k27me3 = (
            h3k27me3
            .drop_duplicates(
                subset="region_id",
                keep="first"
            )
        )

    print(
        f"Unique H3K27me3 regions: {len(h3k27me3)}"
    )


    print("\nJoining H3K27me3 features...")

    feature_cols = [
        "region_id",
        "H3K27me3_mean",
        "H3K27me3_max",
        "H3K27me3_p90",
        "H3K27me3_sum"
    ]


    merged = guides.merge(
        h3k27me3[feature_cols],
        on="region_id",
        how="left",
        validate="many_to_one"
    )


    matched = (
        merged["H3K27me3_mean"]
        .notna()
        .sum()
    )

    coverage = (
        matched / len(merged) * 100
    )


    print("\nIntegration summary:")
    print(
        f"Final guides: {len(merged)}"
    )
    print(
        f"Guides with H3K27me3: {matched}"
    )
    print(
        f"H3K27me3 coverage: {coverage:.3f}%"
    )


    merged.to_csv(
        OUT,
        index=False
    )


    print("\nSaved:")
    print(OUT)


    with open(LOG, "a") as f:
        f.write(
            "\n\n## F3 H3K27me3 integration\n"
            f"- Input: {INPUT}\n"
            f"- Output: {OUT}\n"
            f"- Final guides: {len(merged)}\n"
            f"- H3K27me3 matched guides: {matched}\n"
            f"- H3K27me3 coverage: {coverage:.3f}%\n"
        )


if __name__ == "__main__":
    main()
