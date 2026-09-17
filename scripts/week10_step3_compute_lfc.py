import pandas as pd
import numpy as np

path = "data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx"
PSEUDOCOUNT = 1

def load_set(sheet, header_row, ann_sheet):
    df = pd.read_excel(path, sheet_name=sheet, header=header_row)
    df.columns = ["sgRNA_sequence", "pDNA", "HT29_A", "HT29_B", "HT29_C", "A375_A", "A375_B", "A375_C"]
    df = df.dropna(subset=["sgRNA_sequence"])
    df = df[df["sgRNA_sequence"].str.match(r'^[ACGT]+$', na=False)]
    for c in ["pDNA", "HT29_A", "HT29_B", "HT29_C", "A375_A", "A375_B", "A375_C"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna()

    ann = pd.read_excel(path, sheet_name=ann_sheet)
    ann.columns = ["sgRNA_sequence", "gene_symbol", "gene_id"]
    df = df.merge(ann, on="sgRNA_sequence", how="left")
    return df

setA = load_set("SetA raw reads", 2, "SetA sgRNA annotations")
setB = load_set("SetB raw reads", 1, "SetB sgRNA annotations")

print(f"SetA guías válidas (numéricas, ACGT-only): {len(setA)}")
print(f"SetB guías válidas: {len(setB)}")

def cpm_normalize(df, cols):
    out = df.copy()
    for c in cols:
        total = out[c].sum()
        out[c + "_cpm"] = (out[c] + PSEUDOCOUNT) / (total + PSEUDOCOUNT * len(out)) * 1e6
    return out

count_cols = ["pDNA", "HT29_A", "HT29_B", "HT29_C", "A375_A", "A375_B", "A375_C"]

results = []
for name, df in [("SetA", setA), ("SetB", setB)]:
    df = cpm_normalize(df, count_cols)
    df["lfc_HT29"] = np.log2(df[["HT29_A_cpm","HT29_B_cpm","HT29_C_cpm"]].mean(axis=1) / df["pDNA_cpm"])
    df["lfc_A375"] = np.log2(df[["A375_A_cpm","A375_B_cpm","A375_C_cpm"]].mean(axis=1) / df["pDNA_cpm"])
    df["set"] = name
    results.append(df[["sgRNA_sequence","gene_symbol","gene_id","set","pDNA","lfc_HT29","lfc_A375"]])

combined = pd.concat(results, ignore_index=True)
combined = combined.drop_duplicates(subset="sgRNA_sequence")

out_path = "data/processed/sanson2018_lfc_scores.csv"
combined.to_csv(out_path, index=False)

print(f"\nTotal guías únicas con LFC calculado: {len(combined)}")
print(f"\nDistribución lfc_HT29:\n{combined['lfc_HT29'].describe()}")
print(f"\nDistribución lfc_A375:\n{combined['lfc_A375'].describe()}")
print(f"\nGuardado en: {out_path}")

# Matching contra matriz Phase 1 (K562, Horlbeck)
phase1 = pd.read_csv("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
phase1_seqs = set(phase1["guide_sequence"].str.upper())
combined_seqs = set(combined["sgRNA_sequence"].str.upper())
overlap = phase1_seqs & combined_seqs

print(f"\n--- MATCHING CONTRA MATRIZ PHASE 1 ---")
print(f"Guías en Phase1 (Horlbeck, K562): {len(phase1_seqs)}")
print(f"Guías en Sanson2018 con LFC: {len(combined_seqs)}")
print(f"Overlap por secuencia exacta: {len(overlap)}")
