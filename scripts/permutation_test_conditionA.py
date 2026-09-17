#!/usr/bin/env python
# scripts/permutation_test_conditionA.py
"""
Paired permutation test (10,000 perms) on Spearman rho:
full (chromatin) vs sequence-only, Condition A, full Sanson set.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PRED = Path("results/conditionA_sequence_vs_full_predictions.csv")
OUT = Path("results/conditionA_permutation_test.json")
N_PERM = 10_000
RNG = np.random.default_rng(42)

def spearman_rho(y, p):
    return spearmanr(y, p).correlation

def main():
    df = pd.read_csv(PRED)
    print(f"Loaded: {df.shape}")

    results = {}
    for cell, ycol in [("HT29", "lfc_HT29"), ("A375", "lfc_A375")]:
        sub = df.dropna(subset=[ycol, "pred_seq_only", "pred_full"]).copy()
        y = sub[ycol].values
        p_seq = sub["pred_seq_only"].values
        p_full = sub["pred_full"].values
        n = len(y)

        rho_seq = spearman_rho(y, p_seq)
        rho_full = spearman_rho(y, p_full)
        delta_obs = abs(rho_full) - abs(rho_seq)  # improvement in |rho|
        # Also signed delta (both negative expected)
        delta_signed = rho_full - rho_seq  # more negative = better for CRISPRi

        # Paired permutation: under H0, labels seq/full are exchangeable per guide
        null_deltas = np.empty(N_PERM)
        for i in range(N_PERM):
            swap = RNG.random(n) < 0.5
            a = np.where(swap, p_full, p_seq)
            b = np.where(swap, p_seq, p_full)
            null_deltas[i] = abs(spearman_rho(y, a)) - abs(spearman_rho(y, b))

        # One-sided: is |rho_full| > |rho_seq|?
        pval = (np.sum(null_deltas >= delta_obs) + 1) / (N_PERM + 1)

        results[cell] = {
            "n": int(n),
            "rho_seq_only": float(rho_seq),
            "rho_full": float(rho_full),
            "delta_abs_rho": float(delta_obs),
            "delta_signed": float(delta_signed),
            "pvalue_onesided_abs": float(pval),
            "n_perm": N_PERM,
        }
        print(f"\n=== {cell} ===")
        print(f"  N = {n:,}")
        print(f"  ρ seq-only = {rho_seq:.4f}")
        print(f"  ρ full     = {rho_full:.4f}")
        print(f"  Δ|ρ|       = {delta_obs:.4f}")
        print(f"  p (10k perm, one-sided |ρ|) = {pval:.6f}")

    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nGuardado: {OUT}")
    print("=== PERMUTATION TEST COMPLETADO ===")

if __name__ == "__main__":
    main()