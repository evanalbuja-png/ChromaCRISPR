import pandas as pd
from scipy.stats import spearmanr, pearsonr

sanson = pd.read_csv("data/processed/sanson2018_lfc_scores.csv")
phase1 = pd.read_csv("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")

sanson["seq_upper"] = sanson["sgRNA_sequence"].str.upper()
phase1["seq_upper"] = phase1["guide_sequence"].str.upper()

merged = phase1.merge(sanson, on="seq_upper", how="inner")
merged = merged.drop_duplicates(subset="seq_upper")

print(f"N guías en overlap (tras dedup): {len(merged)}")

for col in ["lfc_HT29", "lfc_A375"]:
    sub = merged.dropna(subset=[col, "efficacy_score"])
    rho, p_s = spearmanr(sub["efficacy_score"], sub[col])
    r, p_p = pearsonr(sub["efficacy_score"], sub[col])
    print(f"\n{col} vs efficacy_score (Horlbeck, K562) — n={len(sub)}")
    print(f"  Spearman rho = {rho:.4f} (p={p_s:.4g})")
    print(f"  Pearson r    = {r:.4f} (p={p_p:.4g})")

merged.to_csv("data/processed/sanson2018_horlbeck_overlap.csv", index=False)
print("\nGuardado overlap detallado en data/processed/sanson2018_horlbeck_overlap.csv")
