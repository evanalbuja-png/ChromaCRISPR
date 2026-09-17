#!/usr/bin/env python
# scripts/fix3b_rebuild_horlbeck_matrix_dedup.py
"""
Misma lógica que Fix-3 pero deduplicando guide_sequence
en ambos lados antes del merge (keep='first').
"""

from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

OLD = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
H3K4_FIXED = Path("data/interim/features/f3_h3k4me3_features_fixed.csv")
H3K27_FIXED = Path("data/interim/features/f3_h3k27me3_features_fixed.csv")
OUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")

def main():
    df = pd.read_csv(OLD)
    print(f"Consolidada original: {df.shape}")
    print(f"guide_sequence únicos: {df['guide_sequence'].nunique()}")

    h3k4 = pd.read_csv(H3K4_FIXED)
    h3k27 = pd.read_csv(H3K27_FIXED)

    # Deduplicar fixed por guide_sequence
    h3k4 = h3k4.drop_duplicates(subset="guide_sequence", keep="first")
    h3k27 = h3k27.drop_duplicates(subset="guide_sequence", keep="first")
    print(f"H3K4me3 fixed únicos: {len(h3k4)}")
    print(f"H3K27me3 fixed únicos: {len(h3k27)}")

    h3k4_cols = ["H3K4me3_mean", "H3K4me3_max", "H3K4me3_p90", "H3K4me3_sum"]
    h3k27_cols = ["H3K27me3_mean", "H3K27me3_max", "H3K27me3_p90", "H3K27me3_sum"]

    old_vals = df[["guide_sequence", "H3K4me3_mean"]].copy()

    df = df.drop(columns=[c for c in h3k4_cols + h3k27_cols if c in df.columns])

    # También deduplicar la consolidada por si acaso
    n_before = len(df)
    df = df.drop_duplicates(subset="guide_sequence", keep="first")
    print(f"Consolidada tras dedup guide_sequence: {n_before} → {len(df)}")

    df = df.merge(h3k4[["guide_sequence"] + h3k4_cols], on="guide_sequence", how="left")
    df = df.merge(h3k27[["guide_sequence"] + h3k27_cols], on="guide_sequence", how="left")

    print(f"\nShape final: {df.shape}")
    print("Missing:")
    for c in h3k4_cols + h3k27_cols:
        print(f"  {c}: {df[c].isna().sum()}")

    cmp = old_vals.merge(df[["guide_sequence", "H3K4me3_mean"]], on="guide_sequence",
                         suffixes=("_old", "_new"))
    rho, _ = spearmanr(cmp["H3K4me3_mean_old"], cmp["H3K4me3_mean_new"], nan_policy="omit")
    print(f"\nSpearman H3K4me3 OLD vs NEW: {rho:.4f}")

    df.to_csv(OUT, index=False)
    print(f"Guardado: {OUT}")
    print("=== LISTO Fix-3b ===")

if __name__ == "__main__":
    main()