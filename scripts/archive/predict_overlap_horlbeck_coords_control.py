#!/usr/bin/env python
# scripts/predict_overlap_horlbeck_coords_control.py
"""
Control: features de cromatina extraídas en coordenadas HORLBECK
para las 526 guías del overlap → predecir con el modelo → correlacionar vs LFC.
"""

from pathlib import Path
import gzip
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MAT_DIR = Path("data/interim/conditionA/overlap_control/matrices")
TABLE = Path("data/interim/conditionA/overlap_control/overlap_guides_horlbeck_coords.csv")
MODEL = Path("models/xgb_phase1_horlbeck_conditionA.joblib")
FEATURES = Path("models/xgb_phase1_feature_list.json")
HORLBECK_FULL = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
OUT = Path("results/conditionA_overlap_horlbeck_coords_control.csv")

def load_gz(path):
    rids, signals = [], []
    with gzip.open(path, "rt") as f:
        f.readline()  # JSON header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rids.append(parts[3])
            bins = np.array([float(x) if x not in (".", "nan", "") else np.nan for x in parts[6:]])
            signals.append(bins)
    return np.array(rids), np.vstack(signals)

def summarize(signal):
    return (np.nanmean(signal, 1), np.nanmax(signal, 1),
            np.nanpercentile(signal, 90, 1), np.nansum(signal, 1))

def main():
    print("=== 1. Cargando tabla de overlap ===")
    guides = pd.read_csv(TABLE)
    print(f"N = {len(guides)}")

    # --------------------------------------------------
    print("\n=== 2. Features de cromatina (coords Horlbeck) ===")
    chrom_dfs = []
    for mark in ["ATAC", "H3K27ac", "H3K4me3", "H3K27me3"]:
        path = MAT_DIR / f"{mark}_overlap_matrix.gz"
        rids, signal = load_gz(path)
        mean, mx, p90, sm = summarize(signal)
        tmp = pd.DataFrame({
            "region_id": rids,
            f"{mark}_mean": mean,
            f"{mark}_max": mx,
            f"{mark}_p90": p90,
            f"{mark}_sum": sm,
        })
        chrom_dfs.append(tmp)
        print(f"  {mark}: {len(tmp)} filas")

    chrom = chrom_dfs[0]
    for d in chrom_dfs[1:]:
        chrom = chrom.merge(d, on="region_id", how="outer")

    # --------------------------------------------------
    print("\n=== 3. Traer F1 + F5 del training (mismas secuencias) ===")
    h = pd.read_csv(HORLBECK_FULL)
    f1f5_cols = ["guide_sequence", "gc_content", "mfe_rnafold", "guide_length",
                 "g_run_max", "poly_t_count", "nearest_tss_distance", "within_promoter_2kb"]
    f1f5 = h[f1f5_cols].drop_duplicates("guide_sequence")
    mat = guides.merge(f1f5, on="guide_sequence", how="left")
    mat = mat.merge(chrom, on="region_id", how="left")

    # --------------------------------------------------
    print("\n=== 4. Predicción ===")
    model = joblib.load(MODEL)
    with open(FEATURES) as f:
        feats = json.load(f)

    X = mat[feats].copy()
    X = X.fillna(X.median())
    mat["pred_efficacy"] = model.predict(X)

    # --------------------------------------------------
    print("\n=== 5. Correlaciones ===")
    print(f"{'Métrica':<45} {'ρ':>8}")
    print("-" * 55)

    # A) Score experimental crudo vs LFC (referencia Semana 10)
    for cell in ["lfc_HT29", "lfc_A375"]:
        rho, _ = spearmanr(mat["efficacy_score"], mat[cell], nan_policy="omit")
        print(f"efficacy_score (crudo) vs {cell:<22} {rho:8.4f}")

    # B) Modelo con features en coords HORLBECK vs LFC
    for cell in ["lfc_HT29", "lfc_A375"]:
        rho, _ = spearmanr(mat["pred_efficacy"], mat[cell], nan_policy="omit")
        print(f"Modelo (coords Horlbeck) vs {cell:<18} {rho:8.4f}")

    # C) Comparación de H3K4me3 ahora (debería ser ~1.0)
    # Traer H3K4me3 original del training
    h3 = h[["guide_sequence", "H3K4me3_mean"]].drop_duplicates("guide_sequence")
    mat2 = mat.merge(h3, on="guide_sequence", how="left", suffixes=("", "_train"))
    rho_h3, _ = spearmanr(mat2["H3K4me3_mean"], mat2["H3K4me3_mean_train"], nan_policy="omit")
    print(f"\nH3K4me3 (re-extraído coords H) vs training:  {rho_h3:8.4f}")

    mat.to_csv(OUT, index=False)
    print(f"\nResultados guardados: {OUT}")
    print("=== CONTROL COMPLETADO ===")

if __name__ == "__main__":
    main()