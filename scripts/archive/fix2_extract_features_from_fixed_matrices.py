#!/usr/bin/env python
# scripts/fix2_extract_features_from_fixed_matrices.py
"""
Lee las matrices fixed de H3K4me3 y H3K27me3,
usa el NAME de deepTools como region_id/clave,
calcula mean/max/p90/sum y guarda CSVs alineados.
"""

from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd

MAT_DIR = Path("data/interim/features/histones_fixed")
BED = Path("data/interim/features/sgRNA_hg38.bed")
OUT_DIR = Path("data/interim/features")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def load_matrix(path):
    """Devuelve DataFrame con region_id (=name) + bins."""
    region_ids = []
    chroms, starts, ends = [], [], []
    signals = []
    with gzip.open(path, "rt") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # chrom, start, end, name, score, strand, bin1...
            chroms.append(parts[0])
            starts.append(int(parts[1]))
            ends.append(int(parts[2]))
            region_ids.append(parts[3])  # NAME — clave correcta
            bins = np.array([float(x) if x not in (".", "nan", "NaN", "") else np.nan
                             for x in parts[6:]])
            signals.append(bins)
    signal = np.vstack(signals)
    return pd.DataFrame({
        "chromosome": chroms,
        "start": starts,
        "end": ends,
        "region_id": region_ids,
        "mean": np.nanmean(signal, axis=1),
        "max": np.nanmax(signal, axis=1),
        "p90": np.nanpercentile(signal, 90, axis=1),
        "sum": np.nansum(signal, axis=1),
    })

def main():
    # BED canónico (para traer guide_sequence)
    print("=== Cargando BED ===")
    bed = pd.read_csv(
        BED, sep="\t", header=None,
        names=["chromosome", "start", "end", "guide_sequence"]
    )
    bed["region_id"] = (
        bed["chromosome"] + ":" + bed["start"].astype(str) + "-" + bed["end"].astype(str)
    )
    print(f"BED regiones: {len(bed):,}")

    for mark in ["H3K4me3", "H3K27me3"]:
        print(f"\n=== {mark} ===")
        path = MAT_DIR / f"{mark}_matrix_fixed.gz"
        df = load_matrix(path)
        print(f"Filas matriz: {len(df):,}")
        print(f"region_id únicos: {df['region_id'].nunique():,}")

        # Prefijos de columnas
        df = df.rename(columns={
            "mean": f"{mark}_mean",
            "max": f"{mark}_max",
            "p90": f"{mark}_p90",
            "sum": f"{mark}_sum",
        })

        # Unir guide_sequence desde BED por region_id
        # El name de deepTools suele ser chrom:start-end o el 4º campo del BED
        # En sgRNA_hg38.bed el 4º campo es guide_sequence, no region_id.
        # computeMatrix usa el name del BED → en este caso es guide_sequence.
        # Verificamos:
        print("Ejemplos region_id (name de matriz):", df["region_id"].head(3).tolist())
        print("Ejemplos guide_sequence (BED):", bed["guide_sequence"].head(3).tolist())

        # Caso A: name == guide_sequence
        if df["region_id"].iloc[0] in set(bed["guide_sequence"].head(1000)):
            print("→ name de matriz = guide_sequence")
            merged = bed.merge(
                df[["region_id", f"{mark}_mean", f"{mark}_max", f"{mark}_p90", f"{mark}_sum"]],
                left_on="guide_sequence", right_on="region_id",
                how="left", suffixes=("", "_mat")
            )
            merged = merged.drop(columns=["region_id"], errors="ignore")
        else:
            # Caso B: name == chrom:start-end
            print("→ name de matriz = chrom:start-end")
            merged = bed.merge(
                df[["region_id", f"{mark}_mean", f"{mark}_max", f"{mark}_p90", f"{mark}_sum"]],
                on="region_id", how="left"
            )

        out = OUT_DIR / f"f3_{mark.lower()}_features_fixed.csv"
        cols = ["chromosome", "start", "end", "guide_sequence",
                f"{mark}_mean", f"{mark}_max", f"{mark}_p90", f"{mark}_sum"]
        cols = [c for c in cols if c in merged.columns]
        merged[cols].to_csv(out, index=False)

        n_ok = merged[f"{mark}_mean"].notna().sum()
        print(f"Guardado: {out}")
        print(f"Guías con {mark}: {n_ok:,} / {len(merged):,} ({100*n_ok/len(merged):.1f}%)")
        print(f"mean stats: mean={merged[f'{mark}_mean'].mean():.3f}, "
              f"std={merged[f'{mark}_mean'].std():.3f}")

    print("\n=== LISTO Fix-2 ===")

if __name__ == "__main__":
    main()