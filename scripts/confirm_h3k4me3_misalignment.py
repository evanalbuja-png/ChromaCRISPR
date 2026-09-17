#!/usr/bin/env python
# scripts/confirm_h3k4me3_misalignment.py
import gzip, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# 1. Leer matriz original usando el NAME de deepTools (correcto)
rids, means = [], []
with gzip.open("data/interim/features/histones/H3K4me3_matrix.gz", "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip().split("\t")
        rids.append(p[3])                    # name real
        sig = np.array([float(x) for x in p[6:]])
        means.append(np.mean(sig))

correct = pd.DataFrame({"region_id": rids, "H3K4me3_mean_correct": means})

# 2. Lo que se guardó (mal alineado)
stored = pd.read_csv("data/interim/features/f3_h3k4me3_features.csv")

m = stored.merge(correct, on="region_id", how="inner")
print(f"Overlap region_id: {len(m):,}")
rho, _ = spearmanr(m["H3K4me3_mean"], m["H3K4me3_mean_correct"])
print(f"Spearman stored vs correct (mismo region_id): {rho:.4f}")
print("Si ρ ≈ 1 → el valor en la matriz es bueno y solo falló el join.")
print("Si ρ bajo → además hay otro problema.")