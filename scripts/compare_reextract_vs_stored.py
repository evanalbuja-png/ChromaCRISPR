#!/usr/bin/env python
# scripts/compare_reextract_vs_stored.py
"""
Compara mean/max/p90/sum re-extraídos vs los guardados en la matriz consolidada.
"""

from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

MAT_DIR = Path("data/interim/conditionA/reextract_horlbeck/matrices")
TABLE = Path("data/interim/conditionA/reextract_horlbeck/horlbeck_consolidated_for_reextract.csv")

def load_gz(path):
    rids, signals = [], []
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rids.append(parts[3])
            bins = np.array([float(x) if x not in (".", "nan", "") else np.nan
                             for x in parts[6:]])
            signals.append(bins)
    return np.array(rids), np.vstack(signals)

def summarize(signal):
    return {
        "mean": np.nanmean(signal, 1),
        "max":  np.nanmax(signal, 1),
        "p90":  np.nanpercentile(signal, 90, 1),
        "sum":  np.nansum(signal, 1),
    }

def main():
    guides = pd.read_csv(TABLE)
    print(f"Guías: {len(guides):,}\n")

    marks = {
        "ATAC": "ATAC",
        "H3K27ac": "H3K27ac",
        "H3K4me3": "H3K4me3",
        "H3K27me3": "H3K27me3",
    }

    print(f"{'Mark':<12} {'Stat':<6} {'Spearman':>10} {'Pearson':>10} {'N':>8}")
    print("-" * 55)

    results = []
    for mark, prefix in marks.items():
        path = MAT_DIR / f"{mark}_reextract.gz"
        rids, signal = load_gz(path)
        stats = summarize(signal)

        tmp = pd.DataFrame({"region_id": rids, **{f"re_{k}": v for k, v in stats.items()}})
        merged = guides.merge(tmp, on="region_id", how="inner")

        for stat in ["mean", "max", "p90", "sum"]:
            stored_col = f"{prefix}_{stat}"
            re_col = f"re_{stat}"
            if stored_col not in merged.columns:
                continue
            mask = merged[stored_col].notna() & merged[re_col].notna()
            if mask.sum() < 10:
                continue
            rho, _ = spearmanr(merged.loc[mask, stored_col], merged.loc[mask, re_col])
            r, _ = pearsonr(merged.loc[mask, stored_col], merged.loc[mask, re_col])
            print(f"{mark:<12} {stat:<6} {rho:10.4f} {r:10.4f} {mask.sum():8d}")
            results.append({
                "mark": mark, "stat": stat,
                "spearman": rho, "pearson": r, "n": int(mask.sum())
            })

        # Distribuciones mean
        print(f"  → stored {prefix}_mean: mean={merged[f'{prefix}_mean'].mean():.3f}, "
              f"std={merged[f'{prefix}_mean'].std():.3f}")
        print(f"  → re-extracted mean:   mean={merged['re_mean'].mean():.3f}, "
              f"std={merged['re_mean'].std():.3f}\n")

    pd.DataFrame(results).to_csv(
        "results/reextract_vs_stored_correlation.csv", index=False)
    print("Guardado: results/reextract_vs_stored_correlation.csv")
    print("=== FIN R3 ===")

if __name__ == "__main__":
    main()