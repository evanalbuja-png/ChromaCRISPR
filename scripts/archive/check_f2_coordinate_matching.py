import pandas as pd


GUIDE_FILE = "data/processed/sgRNA_unified_FINAL.csv"
F2_FILE = "data/interim/features/f2_atac_features.csv"


def make_region_id(chrom, coord):
    return f"{chrom}:{int(coord)}-{int(coord)+1}"


print("Loading files...")

guides = pd.read_csv(
    GUIDE_FILE,
    low_memory=False
)

f2 = pd.read_csv(
    F2_FILE,
    low_memory=False
)


print("\nCreating guide region_id...")

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

guides["coordinate"] = guides["coordinate"].astype(int)

guides["region_id"] = guides.apply(
    lambda x: make_region_id(
        x["chromosome"],
        x["coordinate"]
    ),
    axis=1
)


print("\nChecking exact matches...")

matched = guides[
    guides["region_id"].isin(
        f2["region_id"]
    )
].copy()

unmatched = guides[
    ~guides["region_id"].isin(
        f2["region_id"]
    )
].copy()


print(f"Matched guides: {len(matched)}")
print(f"Unmatched guides: {len(unmatched)}")


print("\n===== MATCHED EXAMPLES =====")

print(
    matched[
        [
            "guide_sequence",
            "chromosome",
            "coordinate",
            "region_id"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


print("\n===== UNMATCHED EXAMPLES =====")

print(
    unmatched[
        [
            "guide_sequence",
            "chromosome",
            "coordinate",
            "region_id"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


print("\n===== Checking possible offsets for unmatched =====")

offsets = [-2, -1, 0, 1, 2]

sample = unmatched.head(20)

for offset in offsets:

    count = 0

    for _, row in sample.iterrows():

        test_id = make_region_id(
            row["chromosome"],
            row["coordinate"] + offset
        )

        if test_id in set(f2["region_id"]):
            count += 1

    print(
        f"Offset {offset}: {count}/20 matches"
    )
