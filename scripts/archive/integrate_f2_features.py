#!/usr/bin/env python3

import pandas as pd
from pathlib import Path


GUIDE_FILE = Path(
    "data/processed/sgRNA_unified_FINAL.csv"
)

F2_FILE = Path(
    "data/interim/features/f2_atac_features.csv"
)

OUTPUT_FILE = Path(
    "data/interim/features/sgRNA_with_f2.csv"
)


def main():

    print("Loading guide dataset...")
    guides = pd.read_csv(
        GUIDE_FILE,
        low_memory=False
    )

    print("Loading F2 ATAC features...")
    f2 = pd.read_csv(
        F2_FILE,
        low_memory=False
    )

    print("\nFiltering guides without valid coordinates...")

    guides_before = len(guides)

    guides = guides[
        guides["coordinate"].notna()
    ].copy()

    guides["coordinate"] = pd.to_numeric(
        guides["coordinate"],
        errors="coerce"
    )

    guides = guides[
        guides["coordinate"].notna()
    ].copy()

    guides["coordinate"] = (
        guides["coordinate"]
        .astype(int)
    )

    excluded = guides_before - len(guides)

    print(
        f"Guides excluded due to invalid coordinates: {excluded}"
    )

    print("\nCreating region_id for guides...")

    print("\nCreating region_id with 1-based to 0-based correction...")

    guides["start"] = guides["coordinate"] - 1
    guides["end"] = guides["start"] + 1

    guides["region_id"] = (
        guides["chromosome"].astype(str)
        + ":"
        + guides["start"].astype(int).astype(str)
        + "-"
        + guides["end"].astype(int).astype(str)
    )

    print("Guide regions:")
    print(len(guides))

    print("\nF2 regions:")
    print(len(f2))

    # Mantener únicamente columnas necesarias de F2
    f2_features = f2[
        [
            "region_id",
            "ATAC_mean",
            "ATAC_max",
            "ATAC_p90",
            "ATAC_sum"
        ]
    ].copy()

    print("\nChecking duplicated F2 region_id...")
    print(
        "Duplicates:",
        f2_features["region_id"].duplicated().sum()
    )

    f2_features = (
        f2_features
        .drop_duplicates(
            subset="region_id",
            keep="first"
        )
    )

    print(
        "Unique F2 regions:",
        len(f2_features)
    )

    print("\nJoining F2 features...")

    merged = guides.merge(
        f2_features,
        on="region_id",
        how="left",
        validate="many_to_one"
    )

    matched = merged["ATAC_mean"].notna().sum()

    print("\nIntegration summary:")
    print(
        f"Guides total: {len(merged)}"
    )
    print(
        f"Guides with F2: {matched}"
    )
    print(
        f"Coverage: {matched/len(merged)*100:.3f}%"
    )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    merged.to_csv(
        OUTPUT_FILE,
        index=False
    )

    print("\nSaved:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
