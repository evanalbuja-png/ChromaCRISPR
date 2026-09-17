#!/usr/bin/env python3
"""
Script: evaluate_rs2_loco.py
Descripción: Cruce y evaluación LOCO (Leave-One-Chromosome-Out) Spearman entre 
             Rule Set 2 y EEP_percentile.
Entorno: chromacrispr-phase1
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr


def run_loco_evaluation():
    matrix_csv = "data/processed/K562_training_matrix_v1.csv"
    rs2_csv = "results/rs2_horlbeck_scored.csv"
    
    if not os.path.exists(rs2_csv):
        print(f"[BLOQUEANTE] No se encuentra el archivo {rs2_csv}.")
        print("Ejecute primero 'python scripts/benchmarking/run_azimuth_rs2.py' en el entorno 'chromacrispr-bench'.")
        sys.exit(1)
        
    if not os.path.exists(matrix_csv):
        print(f"[ERROR] No se encuentra la matriz principal: {matrix_csv}")
        sys.exit(1)

    print(f"[INFO] Cargando matriz K562 y scores RS2...")
    df_matrix = pd.read_csv(matrix_csv)
    df_rs2 = pd.read_csv(rs2_csv)
    
    # Merge exacto por guide_sequence con el subset de Horlbeck
    merged = pd.merge(df_matrix, df_rs2, on='guide_sequence', how='inner')
    print(f"[INFO] Guías mapeadas exitosamente: {len(merged)} de {len(df_rs2)}")
    
    # Filtrar cromosomas canónicos (chr1-22 + chrX)
    valid_chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    merged = merged[merged['chromosome'].isin(valid_chrs)].copy()
    
    # Correlación Global
    rho_global, p_global = spearmanr(merged['rs2_score'], merged['EEP_percentile'])
    print("\n--- MÉTRICAS GLOBALES ---")
    print(f"N evaluado: {len(merged)}")
    print(f"Spearman ρ Global: {rho_global:.4f} (p-value: {p_global:.4e})")
    
    # Correlación LOCO (por cromosoma)
    chromosomes = sorted(merged['chromosome'].unique())
    loco_rhos = []
    
    print("\n--- CORRELACIÓN SPEARMAN LOCO POR CROMOSOMA ---")
    for chrom in chromosomes:
        sub = merged[merged['chromosome'] == chrom]
        if len(sub) >= 10:
            rho_chrom, _ = spearmanr(sub['rs2_score'], sub['EEP_percentile'])
            loco_rhos.append(rho_chrom)
            print(f"  {chrom:5s}: n={len(sub):4d} | Spearman ρ = {rho_chrom:.4f}")
            
    mean_loco = np.mean(loco_rhos)
    sd_loco = np.std(loco_rhos)
    
    print("\n==========================================")
    print(" RESULTADO FINAL BENCHMARK LOCO RS2:")
    print(f" Folds (cromosomas): {len(loco_rhos)}")
    print(f" LOCO Spearman ρ:   {mean_loco:.3f} ± {sd_loco:.3f}")
    print("==========================================")


if __name__ == "__main__":
    run_loco_evaluation()