import pandas as pd


BASE = "data/interim/features/sgRNA_with_f2_f3_repressive.csv"

F1 = "data/interim/features/f1_sequence_features.csv"
F5_CONTEXT = "data/interim/features/f5_genomic_context_features.csv"
F5_TSS = "data/interim/features/f5_tss_features.csv"

OUTPUT = "data/interim/features/sgRNA_with_f1_f5.csv"


print("="*60)
print("ChromaCRISPR Phase 1")
print("Week 7 - F1 + F5 Feature Integration")
print("="*60)


print("\nLoading base dataset...")
base = pd.read_csv(BASE, low_memory=False)

print("Base rows:", len(base))


base["sgrna_id"] = range(len(base))


print("\nLoading F1...")
f1 = pd.read_csv(F1)

f1["sgrna_id"] = range(len(f1))


print("\nLoading F5 genomic context...")
f5 = pd.read_csv(F5_CONTEXT)

f5["sgrna_id"] = range(len(f5))


print("\nLoading F5 TSS...")
f5_tss = pd.read_csv(F5_TSS)


print("\nMerging F1...")
df = base.merge(
    f1.drop(columns=["region_id"]),
    on="sgrna_id",
    how="left",
    validate="one_to_one"
)


print("After F1:", len(df))


print("\nMerging F5 context...")
df = df.merge(
    f5.drop(columns=["region_id"]),
    on="sgrna_id",
    how="left",
    validate="one_to_one"
)


print("After F5 context:", len(df))


print("\nMerging F5 TSS...")
df = df.merge(
    f5_tss.drop(columns=["region_id"]),
    on="sgrna_id",
    how="left",
    validate="one_to_one"
)


print("After F5 TSS:", len(df))


print("\nValidation")

print("Expected rows:", len(base))
print("Final rows:", len(df))

print("\nMissing values:")
print(
    df[
        [
            "gc_content",
            "mfe_rnafold",
            "nearest_tss_distance",
            "within_promoter_2kb"
        ]
    ].isna().sum()
)


duplicates = df["sgrna_id"].duplicated().sum()

print("\nDuplicated sgrna_id:", duplicates)


if len(df) != len(base):
    raise RuntimeError("Row count changed during integration")

if duplicates != 0:
    raise RuntimeError("Duplicate sgrna_id detected")


print("\nSaving:")
print(OUTPUT)

df.to_csv(
    OUTPUT,
    index=False
)

print("\nIntegration completed successfully.")
