#!/usr/bin/env python3
"""
build_sanson_feature_reference.py (v2 - Corregido)
Pipeline bioinformático ChromaCRISPR Phase 1
"""

import pandas as pd
from pathlib import Path

# ==================== CONFIG ====================
WINDOW_FILE = "logs/sanson_window_mapped.tsv"
TSS_FILE = "logs/sanson_tss_resolution_complete.tsv"
GENOME_FILE = "logs/sanson_genome_unique_resolution.tsv"
POSITION_FILE = "logs/sanson_remaining_position_resolution.tsv"

SGRNA_FILE = "data/interim/sgRNA_unified.csv"
OUTPUT = "data/interim/sanson_features_reference.tsv"

print("🚀 Iniciando construcción de referencia Sanson (v2) - ChromaCRISPR Phase 1")

# ==================== 1. CARGA BASE ====================
sg = pd.read_csv(SGRNA_FILE, low_memory=False)
sg = sg[sg["dataset"] == "Sanson2018"].copy().reset_index(drop=True)
print(f"Guías Sanson: {len(sg):,}")

base = sg[["guide_sequence", "gene_symbol"]].copy()

# ==================== 2. CARGA Y NORMALIZACIÓN ====================
def load_and_normalize(file_path: str, source_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t")
    df = df.copy()
    
    # Renombramientos clave
    rename_dict = {
        "original_gene": "gene_symbol",
        "position": "guide_position",
        "coordinate_window_start": "guide_position",
        "tss": "tss_position",
    }
    df.rename(columns=rename_dict, inplace=True)
    
    # Columnas estándar
    if "coordinate_source" not in df.columns:
        df["coordinate_source"] = source_name
    if "status" not in df.columns:
        df["status"] = "promoter" if source_name == "GENCODE_window" else pd.NA
    if "crispri_region" not in df.columns:
        df["crispri_region"] = True
    
    # Cálculos si faltan
    if "offset" not in df.columns and "guide_position" in df.columns and "tss_position" in df.columns:
        df["offset"] = df["guide_position"] - df["tss_position"]
    
    # Columnas finales
    keep_cols = ["guide_sequence", "gene_symbol", "guide_position", 
                 "tss_position", "offset", "status", "coordinate_source", "crispri_region"]
    return df[[c for c in keep_cols if c in df.columns]].copy()


print("\nCargando fuentes...")
window = load_and_normalize(WINDOW_FILE, "GENCODE_window")
tss = load_and_normalize(TSS_FILE, "TSS_recovery")
genome = load_and_normalize(GENOME_FILE, "genome_unique_mapping")
position = load_and_normalize(POSITION_FILE, "position_resolution")

print(f"  → GENCODE_window: {len(window):,}")
print(f"  → TSS: {len(tss):,}")
print(f"  → Genome unique: {len(genome):,}")
print(f"  → Position: {len(position):,}")

# ==================== 3. MERGE SECUENCIAL (prioridad) ====================
merged = base.copy()

for src_name, src_df in [
    ("GENCODE_window", window),
    ("TSS_recovery", tss),
    ("genome_unique_mapping", genome),
    ("position_resolution", position)
]:
    print(f"Merging {src_name}...")
    merged = merged.merge(
        src_df,
        on=["guide_sequence", "gene_symbol"],
        how="left",
        suffixes=("", f"_{src_name}")
    )

# ==================== 4. COALESCE LIMPIO ====================
print("\n🔄 Coalesce por columna...")

critical_cols = ["guide_position", "tss_position", "offset", "status", 
                 "coordinate_source", "crispri_region"]

for col in critical_cols:
    # Recolectar todas las versiones de la columna
    col_versions = [c for c in merged.columns if c == col or c.startswith(f"{col}_")]
    if col_versions:
        merged[col] = (
            merged[col_versions]
            .bfill(axis=1)
            .iloc[:, 0]
            .infer_objects(copy=False)
        )
# Eliminar columnas temporales
cols_to_drop = [c for c in merged.columns if any(c.endswith(f"_{s}") for s in 
               ["GENCODE_window", "TSS_recovery", "genome_unique_mapping", "position_resolution"])]
merged = merged.drop(columns=cols_to_drop, errors="ignore")

# Reordenar
final_cols = ["guide_sequence", "gene_symbol", "coordinate_source", 
              "guide_position", "tss_position", "offset", "status", "crispri_region"]
merged = merged[[c for c in final_cols if c in merged.columns] + 
                [c for c in merged.columns if c not in final_cols]]

# ==================== 5. REPORTE ====================
print("\n" + "="*70)
print("REPORTE FINAL")
print("="*70)
print(f"Total guías: {len(merged):,}")
print(
    f"Controles: {(merged['gene_symbol']=='CONTROL').sum():,}"
)
print("\ncoordinate_source:")
print(merged["coordinate_source"].value_counts(dropna=False))
print("\nstatus:")
print(merged["status"].value_counts(dropna=False))
missing = merged[
    merged["guide_position"].isna()
    &
    (merged["gene_symbol"]!="CONTROL")
]

print(
    f"\nSin posición genómica (excluyendo controles): {len(missing):,}"
)

# ==================== 6. GUARDADO ====================
Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(OUTPUT, sep="\t", index=False)
print(f"\n✅ Guardado en: {OUTPUT}")
print(f"Shape: {merged.shape}")