#!/usr/bin/env python
# scripts/diagnose_offset_by_strand.py

import pandas as pd
from pathlib import Path

h = pd.read_csv("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
s = pd.read_csv("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
ref_s = pd.read_csv("data/processed/sanson2018_reference_FINAL.csv")
ref_h = pd.read_csv("data/processed/horlbeck2016_reference_FINAL.csv")

# Unir por secuencia y traer strand de ambas referencias
m = (h[["guide_sequence", "coordinate", "chromosome"]]
     .merge(s[["guide_sequence", "coordinate", "chromosome", "strand"]],
            on="guide_sequence", suffixes=("_h", "_s"))
     .merge(ref_s[["guide_sequence", "strand", "offset", "status"]],
            on="guide_sequence", how="left", suffixes=("", "_ref"))
     .merge(ref_h[["guide_sequence", "strand"]],
            on="guide_sequence", how="left", suffixes=("", "_href")))

m["delta"] = m["coordinate_h"] - m["coordinate_s"]

print("=== Delta por strand (Sanson) ===")
print(m.groupby("strand")["delta"].value_counts().unstack(fill_value=0).head(20))

print("\n=== Delta medio por strand ===")
print(m.groupby("strand")["delta"].agg(["count", "mean", "median"]))

print("\n=== Crosstab strand_Sanson vs strand_Horlbeck ===")
if "strand_href" in m.columns:
    print(pd.crosstab(m["strand"], m["strand_href"], margins=True))

print("\n=== Muestra de 8 guías (strand + delta + coords) ===")
cols = ["guide_sequence", "strand", "coordinate_h", "coordinate_s", "delta"]
print(m[cols].head(8).to_string())