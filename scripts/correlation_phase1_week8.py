#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

INPUT = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
REPORT = Path("logs/week8_correlation_report.md")

FEATURES = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "H3K27ac_mean",
    "H3K4me3_mean",
    "H3K27me3_mean",
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "nearest_tss_distance",
    "within_promoter_2kb",
]

EXPECTED_ROWS = 155454


def fmt(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.4f}"


def add_section(lines, title, pairs, pearson, spearman):
    lines.append(f"\n## {title}\n")
    lines.append(
        "| Feature 1 | Feature 2 | Pearson r | Spearman rho |\n"
        "|---|---|---:|---:|\n"
    )

    for a, b in pairs:
        lines.append(
            f"| `{a}` | `{b}` | "
            f"{fmt(pearson.loc[a, b])} | "
            f"{fmt(spearman.loc[a, b])} |\n"
        )


def main():

    print("=" * 70)
    print("ChromaCRISPR Phase 1 - Week 8")
    print("Paso 3 - Análisis de correlaciones")
    print("=" * 70)

    df = pd.read_csv(INPUT, low_memory=False)

    print(f"\nRows: {len(df):,}")
    print(f"Columns: {len(df.columns)}")

    if len(df) != EXPECTED_ROWS:
        raise ValueError(
            f"Se esperaban {EXPECTED_ROWS} filas; "
            f"se encontraron {len(df)}."
        )

    missing_columns = [
        c for c in FEATURES
        if c not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Faltan columnas requeridas: {missing_columns}"
        )

    # Convertir explícitamente a numérico.
    for col in FEATURES:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Pearson y Spearman con pares completos.
    pearson = df[FEATURES].corr(
        method="pearson"
    )

    spearman = df[FEATURES].corr(
        method="spearman"
    )

    # ---------------------------------------------------------
    # Strong correlations
    # ---------------------------------------------------------

    strong = []

    for i, a in enumerate(FEATURES):
        for b in FEATURES[i + 1:]:
            r = pearson.loc[a, b]
            rho = spearman.loc[a, b]

            if pd.notna(r) and abs(r) > 0.5:
                strong.append(
                    {
                        "feature_1": a,
                        "feature_2": b,
                        "pearson": r,
                        "spearman": rho,
                        "abs_pearson": abs(r),
                    }
                )

    strong = sorted(
        strong,
        key=lambda x: x["abs_pearson"],
        reverse=True,
    )

    # ---------------------------------------------------------
    # Specific biological relationships
    # ---------------------------------------------------------

    atac_h3k27ac = [
        ("ATAC_mean", "H3K27ac_mean"),
        ("ATAC_max", "H3K27ac_mean"),
        ("ATAC_p90", "H3K27ac_mean"),
    ]

    h3k4_atac_h3k27ac = [
        ("H3K4me3_mean", "ATAC_mean"),
        ("H3K4me3_mean", "ATAC_max"),
        ("H3K4me3_mean", "ATAC_p90"),
        ("H3K4me3_mean", "H3K27ac_mean"),
    ]

    h3k27me3_active = [
        ("H3K27me3_mean", "ATAC_mean"),
        ("H3K27me3_mean", "ATAC_max"),
        ("H3K27me3_mean", "ATAC_p90"),
        ("H3K27me3_mean", "H3K27ac_mean"),
        ("H3K27me3_mean", "H3K4me3_mean"),
    ]

    tss_chromatin = [
        ("nearest_tss_distance", "ATAC_mean"),
        ("nearest_tss_distance", "ATAC_max"),
        ("nearest_tss_distance", "ATAC_p90"),
        ("nearest_tss_distance", "H3K27ac_mean"),
        ("nearest_tss_distance", "H3K4me3_mean"),
        ("nearest_tss_distance", "H3K27me3_mean"),
    ]

    # ---------------------------------------------------------
    # Report
    # ---------------------------------------------------------

    REPORT.parent.mkdir(parents=True, exist_ok=True)

    lines = []

    lines.append("# ChromaCRISPR Phase 1 — Week 8\n")
    lines.append("## Paso 3 — Análisis de correlaciones\n")

    lines.append("### Dataset\n")
    lines.append(
        f"- Archivo: `{INPUT}`\n"
        f"- Filas: {len(df):,}\n"
        f"- Features analizadas: {len(FEATURES)}\n"
        f"- Correlación Pearson: sí\n"
        f"- Correlación Spearman: sí\n"
        f"- Método: correlación por pares completos (`pairwise complete observations`)\n"
    )

    # ---------------------------------------------------------
    # Full Pearson matrix
    # ---------------------------------------------------------

    lines.append("\n## 1. Matriz de correlación Pearson\n")

    lines.append(
        pearson.to_markdown(
            floatfmt=".4f"
        )
    )
    lines.append("\n")

    # ---------------------------------------------------------
    # Full Spearman matrix
    # ---------------------------------------------------------

    lines.append("\n## 2. Matriz de correlación Spearman\n")

    lines.append(
        spearman.to_markdown(
            floatfmt=".4f"
        )
    )
    lines.append("\n")

    # ---------------------------------------------------------
    # Strong correlations
    # ---------------------------------------------------------

    lines.append("\n## 3. Correlaciones fuertes (|Pearson r| > 0.5)\n")

    if strong:
        lines.append(
            "| Feature 1 | Feature 2 | Pearson r | Spearman rho |\n"
            "|---|---|---:|---:|\n"
        )

        for item in strong:
            lines.append(
                f"| `{item['feature_1']}` | "
                f"`{item['feature_2']}` | "
                f"{item['pearson']:.4f} | "
                f"{item['spearman']:.4f} |\n"
            )
    else:
        lines.append(
            "No se encontraron correlaciones con `|Pearson r| > 0.5`.\n"
        )

    # ---------------------------------------------------------
    # Specific relationships
    # ---------------------------------------------------------

    add_section(
        lines,
        "4. ATAC vs H3K27ac",
        atac_h3k27ac,
        pearson,
        spearman,
    )

    add_section(
        lines,
        "5. H3K4me3 vs ATAC / H3K27ac",
        h3k4_atac_h3k27ac,
        pearson,
        spearman,
    )

    add_section(
        lines,
        "6. H3K27me3 vs features activas",
        h3k27me3_active,
        pearson,
        spearman,
    )

    add_section(
        lines,
        "7. Distancia al TSS vs cromatina",
        tss_chromatin,
        pearson,
        spearman,
    )

    # ---------------------------------------------------------
    # Basic interpretation
    # ---------------------------------------------------------

    lines.append("\n## 8. Interpretación automática básica\n")

    if strong:
        lines.append(
            "Las correlaciones con `|Pearson r| > 0.5` se consideran "
            "fuertes para este análisis descriptivo y se listan arriba.\n"
        )
    else:
        lines.append(
            "No se observaron correlaciones lineales fuertes "
            "según el umbral definido (`|Pearson r| > 0.5`).\n"
        )

    # Compare Pearson vs Spearman for strong pairs.
    if strong:
        lines.append(
            "La comparación Pearson/Spearman permite detectar si las "
            "relaciones observadas son principalmente lineales o si "
            "también existe asociación monotónica.\n"
        )

    lines.append(
        "\n> Este análisis es exploratorio. Las correlaciones no implican "
        "causalidad y no se realizan eliminación de features, "
        "selección de variables ni transformación de datos en este paso.\n"
    )

    REPORT.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # Console output
    # ---------------------------------------------------------

    print("\n=== CORRELACIONES FUERTES |Pearson r| > 0.5 ===")

    if strong:
        for item in strong:
            print(
                f"{item['feature_1']} <-> {item['feature_2']}: "
                f"Pearson={item['pearson']:.4f}, "
                f"Spearman={item['spearman']:.4f}"
            )
    else:
        print("Ninguna.")

    print("\n=== ATAC vs H3K27ac ===")
    for a, b in atac_h3k27ac:
        print(
            f"{a} <-> {b}: "
            f"Pearson={pearson.loc[a,b]:.4f}, "
            f"Spearman={spearman.loc[a,b]:.4f}"
        )

    print("\n=== H3K4me3 vs ATAC / H3K27ac ===")
    for a, b in h3k4_atac_h3k27ac:
        print(
            f"{a} <-> {b}: "
            f"Pearson={pearson.loc[a,b]:.4f}, "
            f"Spearman={spearman.loc[a,b]:.4f}"
        )

    print("\n=== H3K27me3 vs features activas ===")
    for a, b in h3k27me3_active:
        print(
            f"{a} <-> {b}: "
            f"Pearson={pearson.loc[a,b]:.4f}, "
            f"Spearman={spearman.loc[a,b]:.4f}"
        )

    print("\n=== TSS distance vs cromatina ===")
    for a, b in tss_chromatin:
        print(
            f"{a} <-> {b}: "
            f"Pearson={pearson.loc[a,b]:.4f}, "
            f"Spearman={spearman.loc[a,b]:.4f}"
        )

    print(f"\nReporte generado: {REPORT}")


if __name__ == "__main__":
    main()
