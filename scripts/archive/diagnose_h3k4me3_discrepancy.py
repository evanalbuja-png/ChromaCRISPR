#!/usr/bin/env python
# scripts/diagnose_h3k4me3_discrepancy.py
"""
Por qué H3K4me3_mean tiene ρ ≈ 0 entre la matriz Horlbeck y la matriz Sanson
para las mismas 526 guías.
"""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr

HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
SANSON = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")

def main():
    h = pd.read_csv(HORLBECK)
    s = pd.read_csv(SANSON)

    m = h.merge(s, on="guide_sequence", how="inner", suffixes=("_h", "_s"))
    print(f"Overlap: {len(m):,} guías\n")

    # 1. Distribuciones
    print("=== Distribuciones H3K4me3_mean ===")
    for label, col in [("Horlbeck (training)", "H3K4me3_mean_h"),
                       ("Sanson (Condition A)", "H3K4me3_mean_s")]:
        x = m[col].dropna()
        print(f"{label}:")
        print(f"  n={len(x)}, mean={x.mean():.3f}, std={x.std():.3f}, "
              f"min={x.min():.3f}, max={x.max():.3f}, "
              f"median={x.median():.3f}")
        print(f"  % ceros: {(x == 0).mean():.1%}")
        print(f"  % NaN original: {m[col].isna().mean():.1%}")

    # 2. Correlación y scatter stats
    rho, _ = spearmanr(m["H3K4me3_mean_h"], m["H3K4me3_mean_s"], nan_policy="omit")
    r, _ = pearsonr(m["H3K4me3_mean_h"].fillna(0), m["H3K4me3_mean_s"].fillna(0))
    print(f"\nSpearman H3K4me3_h vs _s: {rho:.4f}")
    print(f"Pearson  H3K4me3_h vs _s: {r:.4f}")

    # 3. Comparar con ATAC (control)
    print("\n=== Control: ATAC_mean (debería ser ~1.0) ===")
    rho_a, _ = spearmanr(m["ATAC_mean_h"], m["ATAC_mean_s"], nan_policy="omit")
    print(f"Spearman ATAC_h vs _s: {rho_a:.4f}")

    # 4. ¿Las coordenadas coinciden?
    print("\n=== Coordenadas del mismo guide ===")
    if "coordinate_h" in m.columns and "coordinate_s" in m.columns:
        coord_diff = (m["coordinate_h"] - m["coordinate_s"]).abs()
        print(f"Diff coordenada (abs): mean={coord_diff.mean():.1f}, "
              f"max={coord_diff.max()}, % idénticas={(coord_diff==0).mean():.1%}")
    elif "coordinate" in m.columns:
        print("Solo una columna 'coordinate' después del merge — revisar suffixes")

    # 5. Muestra de 10 guías con mayor discrepancia
    m["diff_h3k4"] = (m["H3K4me3_mean_h"] - m["H3K4me3_mean_s"]).abs()
    print("\n=== Top 10 mayores discrepancias H3K4me3 ===")
    cols = ["guide_sequence", "H3K4me3_mean_h", "H3K4me3_mean_s", "diff_h3k4"]
    if "coordinate_h" in m.columns:
        cols += ["coordinate_h", "coordinate_s"]
    print(m.nlargest(10, "diff_h3k4")[cols].to_string())

    # 6. ¿Cuántos H3K4me3 son exactamente 0 en cada set?
    print("\n=== Ceros y valores bajos ===")
    print(f"Horlbeck H3K4me3 == 0: {(m['H3K4me3_mean_h']==0).sum()}")
    print(f"Sanson   H3K4me3 == 0: {(m['H3K4me3_mean_s']==0).sum()}")
    print(f"Horlbeck H3K4me3 < 0.1: {(m['H3K4me3_mean_h']<0.1).sum()}")
    print(f"Sanson   H3K4me3 < 0.1: {(m['H3K4me3_mean_s']<0.1).sum()}")

if __name__ == "__main__":
    main()