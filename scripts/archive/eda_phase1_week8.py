#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

INPUT = Path("data/processed/sgRNA_feature_matrix_phase1.csv")
REPORT = Path("logs/week8_eda_report.md")
PLOTS = Path("logs/week8_eda_plots")

EXPECTED_ROWS = 155454

KEY_FEATURES = [
    "ATAC_mean",
    "H3K27ac_mean",
    "H3K4me3_mean",
    "H3K27me3_mean",
    "gc_content",
    "mfe_rnafold",
    "nearest_tss_distance",
]

F1 = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count",
]

F2 = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum",
]

F3 = [
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

F5 = [
    "nearest_tss_distance",
    "within_promoter_2kb",
]


def check_columns(df, columns, group):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"Faltan columnas de {group}: {missing}"
        )


def fmt(value):
    if pd.isna(value):
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main():

    print("=" * 70)
    print("ChromaCRISPR Phase 1 - Week 8")
    print("Paso 2 - EDA básico")
    print("=" * 70)

    # ---------------------------------------------------------
    # Load matrix
    # ---------------------------------------------------------

    print("\nLoading final Phase 1 matrix...")
    df = pd.read_csv(INPUT, low_memory=False)

    print("Rows:", len(df))
    print("Columns:", len(df.columns))

    if len(df) != EXPECTED_ROWS:
        raise ValueError(
            f"Se esperaban {EXPECTED_ROWS} filas, pero se encontraron {len(df)}."
        )

    check_columns(df, F1, "F1")
    check_columns(df, F2, "F2")
    check_columns(df, F3, "F3")
    check_columns(df, F5, "F5")
    check_columns(df, KEY_FEATURES, "features clave")

    # ---------------------------------------------------------
    # Numeric columns
    # ---------------------------------------------------------

    feature_columns = F1 + F2 + F3 + ["nearest_tss_distance"]

    for col in feature_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # ---------------------------------------------------------
    # Create report
    # ---------------------------------------------------------

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)

    lines = []

    lines.append("# ChromaCRISPR Phase 1 — Week 8 EDA\n")
    lines.append("## Paso 2 — Análisis exploratorio básico\n")

    lines.append("### Dataset\n")
    lines.append(
        f"- Archivo: `{INPUT}`\n"
        f"- Filas: {len(df):,}\n"
        f"- Columnas: {len(df.columns)}\n"
        f"- Features F1: {len(F1)}\n"
        f"- Features F2: {len(F2)}\n"
        f"- Features F3: {len(F3)}\n"
        f"- Features F5: {len(F5)}\n"
    )

    # ---------------------------------------------------------
    # 1. Statistical summary
    # ---------------------------------------------------------

    lines.append("## 1. Resumen estadístico de features\n")

    summary = df[feature_columns].describe(
        percentiles=[0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]
    ).T

    summary["missing"] = df[feature_columns].isna().sum()
    summary["missing_pct"] = (
        summary["missing"] / len(df) * 100
    )

    lines.append(
        "| Feature | Count | Missing | Missing % | Mean | Std | Min | "
        "1% | 5% | 25% | 50% | 75% | 95% | 99% | Max |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
    )

    for col, row in summary.iterrows():
        lines.append(
            f"| `{col}` | {fmt(row['count'])} | "
            f"{int(row['missing'])} | {row['missing_pct']:.3f}% | "
            f"{fmt(row['mean'])} | {fmt(row['std'])} | "
            f"{fmt(row['min'])} | {fmt(row['1%'])} | "
            f"{fmt(row['5%'])} | {fmt(row['25%'])} | "
            f"{fmt(row['50%'])} | {fmt(row['75%'])} | "
            f"{fmt(row['95%'])} | {fmt(row['99%'])} | "
            f"{fmt(row['max'])} |\n"
        )

    # ---------------------------------------------------------
    # 2. Key feature distributions
    # ---------------------------------------------------------

    lines.append("\n## 2. Distribución de features clave\n")

    for col in KEY_FEATURES:
        s = df[col].dropna()

        lines.append(f"### `{col}`\n")

        if len(s) == 0:
            lines.append("- Sin valores disponibles.\n")
            continue

        lines.append(
            f"- N válido: {len(s):,}\n"
            f"- Faltantes: {df[col].isna().sum():,}\n"
            f"- Media: {s.mean():.6g}\n"
            f"- Mediana: {s.median():.6g}\n"
            f"- Desviación estándar: {s.std():.6g}\n"
            f"- Mínimo: {s.min():.6g}\n"
            f"- Máximo: {s.max():.6g}\n"
            f"- P01: {s.quantile(0.01):.6g}\n"
            f"- P05: {s.quantile(0.05):.6g}\n"
            f"- P95: {s.quantile(0.95):.6g}\n"
            f"- P99: {s.quantile(0.99):.6g}\n"
        )

    # ---------------------------------------------------------
    # 3. Promoter percentage
    # ---------------------------------------------------------

    lines.append("\n## 3. Guías dentro del promotor ±2 kb\n")

    promoter = df["within_promoter_2kb"]

    valid_promoter = promoter.dropna()

    promoter_true = (
        valid_promoter.astype(bool).sum()
    )

    promoter_false = (
        (~valid_promoter.astype(bool)).sum()
    )

    promoter_pct = (
        promoter_true / len(valid_promoter) * 100
        if len(valid_promoter) > 0
        else float("nan")
    )

    lines.append(
        f"- Valores válidos: {len(valid_promoter):,}\n"
        f"- `True`: {promoter_true:,}\n"
        f"- `False`: {promoter_false:,}\n"
        f"- Faltantes: {promoter.isna().sum():,}\n"
        f"- **Porcentaje dentro del promotor ±2 kb: {promoter_pct:.3f}%**\n"
    )

    # ---------------------------------------------------------
    # 4. Dataset distribution
    # ---------------------------------------------------------

    lines.append("\n## 4. Distribución por dataset\n")

    dataset_counts = (
        df["dataset"]
        .value_counts(dropna=False)
        .rename_axis("dataset")
        .reset_index(name="n")
    )

    dataset_counts["percentage"] = (
        dataset_counts["n"] / len(df) * 100
    )

    lines.append(
        "| Dataset | N | Porcentaje |\n"
        "|---|---:|---:|\n"
    )

    for _, row in dataset_counts.iterrows():
        dataset_name = (
            "NA" if pd.isna(row["dataset"])
            else str(row["dataset"])
        )

        lines.append(
            f"| `{dataset_name}` | {int(row['n']):,} | "
            f"{row['percentage']:.3f}% |\n"
        )

    # ---------------------------------------------------------
    # Expected datasets check
    # ---------------------------------------------------------

    expected_datasets = {
        "sanson2018",
        "Horlbeck2016",
        "Gasperini2019",
        "replogle2022",
    }

    observed_datasets = set(
        df["dataset"].dropna().astype(str).unique()
    )

    unexpected = sorted(observed_datasets - expected_datasets)
    missing_expected = sorted(expected_datasets - observed_datasets)

    lines.append("\n### Validación de datasets esperados\n")

    if not unexpected and not missing_expected:
        lines.append(
            "- Los cuatro datasets esperados están presentes.\n"
        )
    else:
        if unexpected:
            lines.append(
                f"- Datasets no esperados encontrados: `{unexpected}`\n"
            )
        if missing_expected:
            lines.append(
                f"- Datasets esperados ausentes: `{missing_expected}`\n"
            )

    # ---------------------------------------------------------
    # 5. Missingness summary
    # ---------------------------------------------------------

    lines.append("\n## 5. Valores faltantes\n")

    missing = df.isna().sum()
    missing_pct = missing / len(df) * 100

    lines.append(
        "| Columna | Faltantes | Porcentaje |\n"
        "|---|---:|---:|\n"
    )

    for col in df.columns:
        if missing[col] > 0:
            lines.append(
                f"| `{col}` | {missing[col]:,} | "
                f"{missing_pct[col]:.3f}% |\n"
            )

    # ---------------------------------------------------------
    # 6. Basic observations
    # ---------------------------------------------------------

    lines.append("\n## 6. Observaciones automáticas\n")

    for col in KEY_FEATURES:
        s = df[col].dropna()

        if len(s) == 0:
            continue

        skew = s.skew()

        lines.append(
            f"- `{col}`: skewness = {skew:.4f}; "
            f"mediana = {s.median():.6g}; "
            f"media = {s.mean():.6g}.\n"
        )

    lines.append(
        "\n> Este EDA es descriptivo. No se realizan transformaciones, "
        "imputaciones ni eliminación de outliers en este paso.\n"
    )

    REPORT.write_text(
        "".join(lines),
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # Console summary
    # ---------------------------------------------------------

    print("\n=== DATASET ===")
    print(f"Rows: {len(df):,}")
    print(f"Columns: {len(df.columns)}")

    print("\n=== KEY FEATURES ===")
    print(
        df[KEY_FEATURES].describe(
            percentiles=[0.05, 0.50, 0.95]
        ).T.to_string()
    )

    print("\n=== WITHIN PROMOTER ±2 KB ===")
    print(f"True: {promoter_true:,}")
    print(f"False: {promoter_false:,}")
    print(f"Missing: {promoter.isna().sum():,}")
    print(f"Percentage True: {promoter_pct:.3f}%")

    print("\n=== DATASETS ===")
    print(dataset_counts.to_string(index=False))

    print("\n=== REPORT ===")
    print(REPORT)


if __name__ == "__main__":
    main()
