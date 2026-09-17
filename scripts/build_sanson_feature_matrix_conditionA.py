#!/usr/bin/env python
# scripts/build_sanson_feature_matrix_conditionA_v2.py
"""
Versión corregida: parsea los .gz de computeMatrix (formato correcto).
"""

from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd
import subprocess
import tempfile
import os

# --------------------------------------------------
GUIDES = Path("data/interim/conditionA/sanson_guides_with_coords.csv")
MAT_DIR = Path("data/interim/conditionA/matrices")
OUT_MATRIX = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
TSS_BED = Path("data/interim/f5_tss.bed")

FEATURES_ORDER = [
    "gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count",
    "ATAC_mean", "ATAC_max", "ATAC_p90", "ATAC_sum",
    "H3K27ac_mean", "H3K27ac_max", "H3K27ac_p90", "H3K27ac_sum",
    "H3K4me3_mean", "H3K4me3_max", "H3K4me3_p90", "H3K4me3_sum",
    "H3K27me3_mean", "H3K27me3_max", "H3K27me3_p90", "H3K27me3_sum",
    "nearest_tss_distance", "within_promoter_2kb",
]

def load_deeptools_gz(path):
    """Lee el .gz de computeMatrix. Devuelve region_id y matriz de bins."""
    region_ids = []
    signals = []
    with gzip.open(path, "rt") as f:
        # Primera línea = JSON metadata
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # chrom, start, end, name, score, strand, bin1, bin2, ...
            rid = parts[3]
            bins = np.array([float(x) if x not in (".", "nan", "NaN", "") else np.nan
                             for x in parts[6:]])
            region_ids.append(rid)
            signals.append(bins)
    signal = np.vstack(signals)
    return np.array(region_ids), signal

def summarize(signal):
    mean = np.nanmean(signal, axis=1)
    mx   = np.nanmax(signal, axis=1)
    p90  = np.nanpercentile(signal, 90, axis=1)
    sm   = np.nansum(signal, axis=1)
    return mean, mx, p90, sm

def compute_mfe_batch(sequences):
    """Calcula MFE con RNAfold (CLI) de forma eficiente."""
    unique = list(dict.fromkeys(sequences))  # preserve order, unique
    print(f"  Secuencias únicas para MFE: {len(unique):,}")

    # Escribir fasta temporal
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        for i, seq in enumerate(unique):
            fh.write(f">s{i}\n{seq.upper()}\n")
        fa_path = fh.name

    # Llamar RNAfold
    out_path = fa_path + ".out"
    cmd = f"RNAfold --noPS < {fa_path} > {out_path}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    # Parsear salida
    mfe_map = {}
    with open(out_path) as f:
        lines = f.readlines()
    for i in range(0, len(lines), 3):
        # línea 0: >sN
        # línea 1: secuencia
        # línea 2: estructura (mfe)
        header = lines[i].strip()
        mfe_line = lines[i+2].strip()
        # formato: ..... ( -3.20)
        mfe = float(mfe_line.split()[-1].strip("()"))
        idx = int(header[2:])  # sN → N
        mfe_map[unique[idx]] = mfe

    os.unlink(fa_path)
    os.unlink(out_path)
    return mfe_map

def main():
    print("=== 1. Cargando guías ===")
    guides = pd.read_csv(GUIDES)
    print(f"Guías: {len(guides):,}")

    # --------------------------------------------------
    print("\n=== 2. Matrices de cromatina (.gz) ===")
    matrices = {
        "ATAC":     MAT_DIR / "atac_sanson_matrix.gz",
        "H3K27ac":  MAT_DIR / "h3k27ac_sanson_matrix.gz",
        "H3K4me3":  MAT_DIR / "h3k4me3_sanson_matrix.gz",
        "H3K27me3": MAT_DIR / "h3k27me3_sanson_matrix.gz",
    }

    chrom_list = []
    for name, path in matrices.items():
        print(f"  → {name} ({path.name})")
        rids, signal = load_deeptools_gz(path)
        mean, mx, p90, sm = summarize(signal)
        tmp = pd.DataFrame({
            "region_id": rids,
            f"{name}_mean": mean,
            f"{name}_max": mx,
            f"{name}_p90": p90,
            f"{name}_sum": sm,
        })
        chrom_list.append(tmp)
        print(f"     filas: {len(tmp):,}")

    chrom = chrom_list[0]
    for df in chrom_list[1:]:
        chrom = chrom.merge(df, on="region_id", how="outer")
    print(f"Cromatina final: {chrom.shape}")

    # --------------------------------------------------
    print("\n=== 3. F1 – Sequence features + MFE ===")
    f1 = guides[["region_id", "guide_sequence"]].copy()
    seqs = f1["guide_sequence"].astype(str).str.upper()
    f1["gc_content"] = (seqs.str.count("G") + seqs.str.count("C")) / seqs.str.len()
    f1["guide_length"] = seqs.str.len()
    f1["poly_t_count"] = seqs.str.count("TTTT")

    def max_g_run(s):
        mx = cur = 0
        for b in s:
            if b == "G":
                cur += 1
                mx = max(mx, cur)
            else:
                cur = 0
        return mx
    f1["g_run_max"] = seqs.apply(max_g_run)

    # MFE via RNAfold CLI
    try:
        mfe_map = compute_mfe_batch(seqs.tolist())
        f1["mfe_rnafold"] = f1["guide_sequence"].map(mfe_map)
        print(f"  MFE calculado. Cobertura: {f1['mfe_rnafold'].notna().mean():.1%}")
    except Exception as e:
        print(f"  ERROR calculando MFE: {e}")
        print("  → mfe_rnafold = 0.0 (temporal)")
        f1["mfe_rnafold"] = 0.0

    # --------------------------------------------------
    print("\n=== 4. F5 – TSS distance ===")
    if TSS_BED.exists():
        tss = pd.read_csv(TSS_BED, sep="\t", header=None,
                          names=["chrom", "start", "end", "name", "score", "strand"])
        print(f"  TSS cargados: {len(tss):,}")
        f5_parts = []
        for chrom_name, grp in guides.groupby("chromosome"):
            tss_chr = tss[tss["chrom"] == chrom_name]
            if len(tss_chr) == 0:
                tmp = grp[["region_id"]].copy()
                tmp["nearest_tss_distance"] = np.nan
                tmp["within_promoter_2kb"] = 0
            else:
                tss_pos = tss_chr["start"].values.astype(float)
                coords = grp["coordinate"].values.astype(float)
                dists = np.min(np.abs(coords[:, None] - tss_pos[None, :]), axis=1)
                tmp = grp[["region_id"]].copy()
                tmp["nearest_tss_distance"] = dists
                tmp["within_promoter_2kb"] = (dists <= 2000).astype(int)
            f5_parts.append(tmp)
        f5 = pd.concat(f5_parts, ignore_index=True)
    else:
        print("  WARNING: no se encontró f5_tss.bed")
        f5 = guides[["region_id"]].copy()
        f5["nearest_tss_distance"] = np.nan
        f5["within_promoter_2kb"] = 0

    # --------------------------------------------------
    print("\n=== 5. Merge final ===")
    mat = guides.merge(chrom, on="region_id", how="left")
    mat = mat.merge(f1.drop(columns=["guide_sequence"]), on="region_id", how="left")
    mat = mat.merge(f5, on="region_id", how="left")

    meta = ["region_id", "guide_sequence", "gene_symbol", "gene_id",
            "chromosome", "coordinate", "strand", "lfc_HT29", "lfc_A375"]
    cols = [c for c in meta + FEATURES_ORDER if c in mat.columns]
    mat = mat[cols]

    mat.to_csv(OUT_MATRIX, index=False)
    print(f"\nMatriz escrita: {OUT_MATRIX}")
    print(f"Shape: {mat.shape}")
    print("\nMissing por feature (top):")
    print(mat[FEATURES_ORDER].isna().sum().sort_values(ascending=False).head(12))
    print("\n=== LISTO ===")

if __name__ == "__main__":
    main()