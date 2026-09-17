#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

BASE = Path("data/interim/features/sgRNA_with_f2_f3_repressive.csv")
F1 = Path("data/interim/features/f1_sequence_features.csv")
F5_CONTEXT = Path("data/interim/features/f5_genomic_context_features.csv")
F5_TSS = Path("data/interim/features/f5_tss_features.csv")

OUT = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
LOG = Path("logs/week8_log.md")

EXPECTED_ROWS = 155454

ID_COLS = [
    "sgrna_id",
    "dataset",
    "experiment",
    "guide_sequence",
    "gene_symbol",
    "gene_id",
    "chromosome",
    "coordinate",
    "strand",
    "region_id",
]

F1_COLS = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count",
]

F2_COLS = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum",
]

F3_COLS = [
    "H3K27ac_mean",
    "H3K27ac_max",
    "H3K27ac_p90",
    "H3K27ac_sum",
    "H3K4me3_mean",
    "H3K4me3_max",
    "H3K4me3_p90",
    "H3K4me3_sum",
    "H3K27me3_mean",
    "H3K27me3_max",
    "H3K27me3_p90",
    "H3K27me3_sum",
]

F5_COLS = [
    "nearest_tss_gene",
    "nearest_tss_distance",
    "within_promoter_2kb",
]

ALL_COLS = ID_COLS + F1_COLS + F2_COLS + F3_COLS + F5_COLS


def check_columns(df, required, name):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name}: faltan columnas requeridas: {missing}"
        )


def main():

    print("=" * 60)
    print("ChromaCRISPR Phase 1")
    print("Week 8 - Final Phase 1 Feature Matrix")
    print("=" * 60)

    # ---------------------------------------------------------
    # Load datasets
    # ---------------------------------------------------------

    print("\nLoading base dataset...")
    base = pd.read_csv(BASE, low_memory=False)
    print("Base rows:", len(base))

    print("\nLoading F1...")
    f1 = pd.read_csv(F1, low_memory=False)
    print("F1 rows:", len(f1))

    print("\nLoading F5 genomic context...")
    f5_context = pd.read_csv(F5_CONTEXT, low_memory=False)
    print("F5 context rows:", len(f5_context))

    print("\nLoading F5 TSS...")
    f5_tss = pd.read_csv(F5_TSS, low_memory=False)
    print("F5 TSS rows:", len(f5_tss))

    # ---------------------------------------------------------
    # Validate row counts
    # ---------------------------------------------------------

    for name, df in [
        ("BASE", base),
        ("F1", f1),
        ("F5_CONTEXT", f5_context),
        ("F5_TSS", f5_tss),
    ]:
        if len(df) != EXPECTED_ROWS:
            raise ValueError(
                f"{name}: {len(df)} filas; se esperaban {EXPECTED_ROWS}."
            )

    # ---------------------------------------------------------
    # Reconstruct sgrna_id exactly as Week 7 integration did
    # ---------------------------------------------------------

    print("\nAssigning positional sgrna_id...")

    base["sgrna_id"] = range(len(base))
    f1["sgrna_id"] = range(len(f1))
    f5_context["sgrna_id"] = range(len(f5_context))

    # F5 TSS already contains sgrna_id.
    if "sgrna_id" not in f5_tss.columns:
        raise ValueError(
            "F5 TSS no contiene sgrna_id; no se modificará automáticamente."
        )

    # ---------------------------------------------------------
    # Validate IDs
    # ---------------------------------------------------------

    for name, df in [
        ("BASE", base),
        ("F1", f1),
        ("F5_CONTEXT", f5_context),
        ("F5_TSS", f5_tss),
    ]:
        duplicated = df["sgrna_id"].duplicated().sum()

        if duplicated:
            raise ValueError(
                f"{name}: {duplicated} sgrna_id duplicados."
            )

        if df["sgrna_id"].nunique() != EXPECTED_ROWS:
            raise ValueError(
                f"{name}: sgrna_id no contiene {EXPECTED_ROWS} valores únicos."
            )

    # ---------------------------------------------------------
    # Check required columns
    # ---------------------------------------------------------

    check_columns(
        base,
        [
            "sgrna_id",
            "dataset",
            "experiment",
            "guide_sequence",
            "gene_symbol",
            "gene_id",
            "chromosome",
            "coordinate",
            "strand",
            "region_id",
        ] + F2_COLS + F3_COLS,
        "BASE",
    )

    check_columns(
        f1,
        ["sgrna_id"] + F1_COLS,
        "F1",
    )

    check_columns(
        f5_context,
        ["sgrna_id"],
        "F5_CONTEXT",
    )

    check_columns(
        f5_tss,
        ["sgrna_id"] + F5_COLS,
        "F5_TSS",
    )

    # ---------------------------------------------------------
    # Select only required columns
    # ---------------------------------------------------------

    base_sel = base[
        [
            "sgrna_id",
            "dataset",
            "experiment",
            "guide_sequence",
            "gene_symbol",
            "gene_id",
            "chromosome",
            "coordinate",
            "strand",
            "region_id",
        ]
        + F2_COLS
        + F3_COLS
    ].copy()

    f1_sel = f1[
        ["sgrna_id"] + F1_COLS
    ].copy()

    f5_tss_sel = f5_tss[
        ["sgrna_id"] + F5_COLS
    ].copy()

    # ---------------------------------------------------------
    # Merge F1
    # ---------------------------------------------------------

    print("\nMerging F1...")

    matrix = base_sel.merge(
        f1_sel,
        on="sgrna_id",
        how="left",
        validate="one_to_one",
    )

    print("After F1:", len(matrix))

    # ---------------------------------------------------------
    # Merge F5 TSS
    # ---------------------------------------------------------

    print("\nMerging F5 TSS...")

    matrix = matrix.merge(
        f5_tss_sel,
        on="sgrna_id",
        how="left",
        validate="one_to_one",
    )

    print("After F5 TSS:", len(matrix))

    # ---------------------------------------------------------
    # Final column order
    # ---------------------------------------------------------

    matrix = matrix[ALL_COLS]

    # ---------------------------------------------------------
    # Final validation
    # ---------------------------------------------------------

    duplicates = matrix["sgrna_id"].duplicated().sum()

    print("\n" + "=" * 60)
    print("FINAL VALIDATION")
    print("=" * 60)

    print("Rows:", len(matrix))
    print("Columns:", len(matrix.columns))
    print("Unique sgrna_id:", matrix["sgrna_id"].nunique())
    print("Duplicated sgrna_id:", duplicates)

    if len(matrix) != EXPECTED_ROWS:
        raise RuntimeError(
            f"Final row count is {len(matrix)}, expected {EXPECTED_ROWS}."
        )

    if duplicates != 0:
        raise RuntimeError(
            "Duplicate sgrna_id detected."
        )

    if matrix["sgrna_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError(
            "sgrna_id is not unique."
        )

    # ---------------------------------------------------------
    # Missing values
    # ---------------------------------------------------------

    print("\nMissing values by column:")
    missing = matrix.isna().sum()

    print(missing.to_string())

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    OUT.parent.mkdir(parents=True, exist_ok=True)

    matrix.to_csv(
        OUT,
        index=False,
    )

    print("\nSaved:")
    print(OUT)

    # ---------------------------------------------------------
    # Update log
    # ---------------------------------------------------------

    LOG.parent.mkdir(parents=True, exist_ok=True)

    with LOG.open("a", encoding="utf-8") as fh:
        fh.write("\n## Paso 1 — Matriz final Phase 1\n\n")
        fh.write(
            "- Construcción de matriz final unificada.\n"
        )
        fh.write(
            f"- Filas: {len(matrix)}\n"
        )
        fh.write(
            f"- Columnas: {len(matrix.columns)}\n"
        )
        fh.write(
            f"- `sgrna_id` únicos: {matrix['sgrna_id'].nunique()}\n"
        )
        fh.write(
            f"- `sgrna_id` duplicados: {duplicates}\n"
        )
        fh.write(
            f"- Archivo: `{OUT}`\n\n"
        )
        fh.write("### Valores faltantes por columna\n\n")

        for col, value in missing.items():
            fh.write(f"- `{col}`: {value}\n")

        fh.write("\n")

    print("\nLog updated:")
    print(LOG)

    print("\nPaso 1 completado.")


if __name__ == "__main__":
    main()
