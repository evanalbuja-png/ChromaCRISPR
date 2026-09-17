#!/usr/bin/env python
# scripts/fix3_rebuild_horlbeck_matrix.py
"""
Reemplaza H3K4me3_* y H3K27me3_* corruptos en la matriz consolidada
por los valores de f3_*_features_fixed.csv (alineados por guide_sequence).
"""

from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

OLD = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
H3K4_FIXED = Path("data/interim/features/f3_h3k4me3_features_fixed.csv")
H3K27_FIXED = Path("data/interim/features/f3_h3k27me3_features_fixed.csv")
OUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")

def main():
    print("=== Cargando matriz consolidada (vieja) ===")
    df = pd.read_csv(OLD)
    print(f"Shape: {df.shape}")
    print(f"Columnas: {df.columns.tolist()}")

    print("\n=== Cargando features fijas ===")
    h3k4 = pd.read_csv(H3K4_FIXED)
    h3k27 = pd.read_csv(H3K27_FIXED)
    print(f"H3K4me3 fixed: {h3k4.shape}")
    print(f"H3K27me3 fixed: {h3k27.shape}")

    # Columnas a reemplazar
    h3k4_cols = ["H3K4me3_mean", "H3K4me3_max", "H3K4me3_p90", "H3K4me3_sum"]
    h3k27_cols = ["H3K27me3_mean", "H3K27me3_max", "H3K27me3_p90", "H3K27me3_sum"]

    # Guardar valores viejos para comparación
    old_h3k4 = df[["guide_sequence"] + h3k4_cols].copy() if all(c in df.columns for c in h3k4_cols) else None

    # Quitar columnas corruptas
    df = df.drop(columns=[c for c in h3k4_cols + h3k27_cols if c in df.columns])

    # Merge por guide_sequence
    df = df.merge(
        h3k4[["guide_sequence"] + h3k4_cols],
        on="guide_sequence",
        how="left"
    )
    df = df.merge(
        h3k27[["guide_sequence"] + h3k27_cols],
        on="guide_sequence",
        how="left"
    )

    print(f"\nShape tras merge: {df.shape}")
    print("Missing tras fix:")
    for c in h3k4_cols + h3k27_cols:
        n = df[c].isna().sum()
        print(f"  {c}: {n} ({100*n/len(df):.1f}%)")

    # Comparación viejo vs nuevo (solo overlap)
    if old_h3k4 is not None:
        cmp = old_h3k4.merge(
            df[["guide_sequence", "H3K4me3_mean"]],
            on="guide_sequence",
            suffixes=("_old", "_new")
        )
        # La columna nueva puede haberse renombrado
        new_col = "H3K4me3_mean_new" if "H3K4me3_mean_new" in cmp.columns else "H3K4me3_mean"
        old_col = "H3K4me3_mean_old" if "H3K4me3_mean_old" in cmp.columns else "H3K4me3_mean"
        if old_col in cmp.columns and new_col in cmp.columns:
            rho, _ = spearmanr(cmp[old_col], cmp[new_col], nan_policy="omit")
            print(f"\nSpearman H3K4me3_mean OLD vs NEW: {rho:.4f}  (esperado ~0.02 si el bug era real)")

    df.to_csv(OUT, index=False)
    print(f"\nMatriz FIXED guardada: {OUT}")
    print(f"Shape final: {df.shape}")
    print("=== LISTO Fix-3 ===")

if __name__ == "__main__":
    main()