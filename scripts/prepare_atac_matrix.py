#!/usr/bin/env python3
"""
prepare_atac_matrix.py

ChromaCRISPR Phase 1 - Semana 5 - Feature F2 (accesibilidad cromatinica, ATAC-seq K562)

Prepara un BED valido a partir de sgRNA_unified.csv (hg38), lo valida
exhaustivamente, y ejecuta deepTools computeMatrix (reference-point, center)
contra el bigWig de ATAC-seq K562 de ENCODE.

Uso:
    conda activate chromacrispr-phase1
    python scripts/prepare_atac_matrix.py

Requisitos: pandas, numpy, pathlib, subprocess (stdlib + pandas/numpy ya
disponibles en el ambiente). No se instalan paquetes nuevos.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Configuracion de rutas (relativas a la raiz del proyecto)
# --------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

INPUT_CSV = PROJECT_ROOT / "data" / "processed" / "sgRNA_unified_FINAL.csv"
OUTPUT_BED = PROJECT_ROOT / "data" / "interim" / "features" / "sgRNA_hg38.bed"
ATAC_BIGWIG = PROJECT_ROOT / "data" / "raw" / "encode_k562" / "ATAC-seq_K562.bigWig"
OUTPUT_MATRIX = PROJECT_ROOT / "data" / "interim" / "features" / "atac_matrix.gz"

REQUIRED_COLUMNS = [
    "guide_sequence",
    "chromosome",
    "coordinate",
    "genome_build",
]

CHR_PATTERN = re.compile(r"^chr([0-9]{1,2}|[XYM]|MT)$")


# --------------------------------------------------------------------------
# Paso 1: carga y validacion del CSV
# --------------------------------------------------------------------------
def load_and_filter(csv_path: Path) -> pd.DataFrame:
    print(f"[1/6] Leyendo dataset unificado: {csv_path}")

    if not csv_path.exists():
        raise FileNotFoundError(f"No se encontro el archivo de entrada: {csv_path}")

    df = pd.read_csv(csv_path)

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Faltan columnas requeridas en el CSV: {missing_cols}. "
            f"Columnas disponibles: {list(df.columns)}"
        )

    n_total = len(df)
    print(f"    Total de filas leidas: {n_total:,}")

    # Filtrar unicamente hg38
    df = df[df["genome_build"] == "hg38"].copy()
    n_hg38 = len(df)
    print(f"    Filas con genome_build == 'hg38': {n_hg38:,}")

    if n_hg38 == 0:
        raise ValueError(
            "Ninguna fila quedo tras filtrar genome_build == 'hg38'. "
            "Verifica el valor exacto de esta columna en el CSV."
        )

    return df


# --------------------------------------------------------------------------
# Paso 2: limpieza y validacion de coordenadas / cromosomas
# --------------------------------------------------------------------------
def clean_and_validate(df: pd.DataFrame) -> pd.DataFrame:
    print("[2/6] Validando y limpiando coordenadas y cromosomas")

    n0 = len(df)

    # chromosome no vacio / no NaN
    df = df[df["chromosome"].notna()].copy()
    df["chromosome"] = df["chromosome"].astype(str).str.strip()
    df = df[df["chromosome"] != ""].copy()

    # Normalizar prefijo chr si falta (defensivo, no oculta errores: se reporta)
    n_missing_prefix = (~df["chromosome"].str.startswith("chr")).sum()
    if n_missing_prefix > 0:
        print(
            f"    AVISO: {n_missing_prefix} filas sin prefijo 'chr' en chromosome; "
            "se antepone 'chr' automaticamente."
        )
        df.loc[~df["chromosome"].str.startswith("chr"), "chromosome"] = (
            "chr" + df.loc[~df["chromosome"].str.startswith("chr"), "chromosome"]
        )

    # coordinate numerico -> coercion; lo no numerico se vuelve NaN y se elimina
    df["coordinate"] = pd.to_numeric(df["coordinate"], errors="coerce")
    df = df[df["coordinate"].notna()].copy()

    # guide_sequence no vacio
    df = df[df["guide_sequence"].notna()].copy()
    df["guide_sequence"] = df["guide_sequence"].astype(str).str.strip()
    df = df[df["guide_sequence"] != ""].copy()

    # Eliminar cualquier NaN residual en columnas clave
    df = df.dropna(subset=["chromosome", "coordinate", "guide_sequence"]).copy()

    # Coordenadas deben ser positivas y convertibles a entero
    df = df[df["coordinate"] >= 0].copy()
    df["coordinate"] = df["coordinate"].astype(np.int64)

    n1 = len(df)
    print(f"    Filas validas tras limpieza: {n1:,} (descartadas: {n0 - n1:,})")

    if n1 == 0:
        raise ValueError("No quedaron filas validas tras la limpieza de coordenadas.")

    # Validar formato de cromosoma (chr1..chr22, chrX, chrY, chrM/chrMT)
    invalid_chr_mask = ~df["chromosome"].apply(lambda c: bool(CHR_PATTERN.match(c)))
    n_invalid_chr = invalid_chr_mask.sum()
    if n_invalid_chr > 0:
        examples = df.loc[invalid_chr_mask, "chromosome"].unique()[:5]
        print(
            f"    AVISO: {n_invalid_chr} filas con cromosoma en formato no estandar "
            f"(ejemplos: {list(examples)}). Se eliminan para evitar errores en deepTools."
        )
        df = df[~invalid_chr_mask].copy()

    n2 = len(df)
    if n2 == 0:
        raise ValueError(
            "No quedaron filas con formato de cromosoma valido (chr1-22, chrX, chrY, chrM)."
        )

    return df


# --------------------------------------------------------------------------
# Paso 3: construccion del BED (0-based, half-open)
# --------------------------------------------------------------------------
def build_bed(df: pd.DataFrame) -> pd.DataFrame:
    print("[3/6] Construyendo coordenadas BED (0-based, half-open)")

    # Se asume que 'coordinate' en el CSV es 1-based (convencion comun en
    # anotaciones de guias). BED usa start 0-based, end 1-based (exclusivo).
    # start = coordinate - 1 ; end = coordinate
    bed = pd.DataFrame()
    bed["chrom"] = df["chromosome"].astype(str)
    bed["start"] = (df["coordinate"] - 1).astype(np.int64)
    bed["end"] = df["coordinate"].astype(np.int64)
    bed["name"] = df["guide_sequence"].astype(str)

    # Nunca debe haber start negativo
    n_negative = (bed["start"] < 0).sum()
    if n_negative > 0:
        print(f"    AVISO: {n_negative} filas con start < 0 tras conversion; se descartan.")
        bed = bed[bed["start"] >= 0].copy()

    if bed.empty:
        raise ValueError("El BED resultante quedo vacio tras construir coordenadas.")

    # Orden estandar por cromosoma y posicion (recomendado por deepTools)
    bed = bed.sort_values(by=["chrom", "start"]).reset_index(drop=True)

    return bed


# --------------------------------------------------------------------------
# Paso 4: escritura del BED con formato estricto
# --------------------------------------------------------------------------
def write_bed(bed: pd.DataFrame, out_path: Path) -> None:
    print(f"[4/6] Escribiendo BED: {out_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    bed.to_csv(
        out_path,
        sep="\t",
        header=False,
        index=False,
        columns=["chrom", "start", "end", "name"],
    )

    print(f"    Escritas {len(bed):,} regiones.")


# --------------------------------------------------------------------------
# Paso 5: validaciones manuales antes de deepTools
# --------------------------------------------------------------------------
def validate_bed_file(bed_path: Path) -> None:
    print(f"[5/6] Validando archivo BED final: {bed_path}")

    if not bed_path.exists() or bed_path.stat().st_size == 0:
        raise RuntimeError(f"El BED generado no existe o esta vacio: {bed_path}")

    with open(bed_path, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    n_lines = len(lines)
    print(f"    Numero de regiones en el BED: {n_lines:,}")
    print("    Primeras 5 lineas:")
    for line in lines[:5]:
        print(f"      {line!r}")

    errors = []
    for i, line in enumerate(lines[:2000]):  # muestreo defensivo para no ser O(n) pesado
        fields = line.split("\t")
        if len(fields) != 4:
            errors.append(f"Linea {i}: numero de columnas != 4 ({len(fields)}): {line!r}")
            continue
        chrom, start, end, name = fields
        if not CHR_PATTERN.match(chrom):
            errors.append(f"Linea {i}: cromosoma con formato invalido: {chrom!r}")
        try:
            start_i, end_i = int(start), int(end)
        except ValueError:
            errors.append(f"Linea {i}: start/end no son enteros: {start!r}, {end!r}")
            continue
        if start_i < 0 or end_i <= start_i:
            errors.append(f"Linea {i}: coordenadas invalidas (start={start_i}, end={end_i})")
        if not name:
            errors.append(f"Linea {i}: nombre de guia vacio")

    if errors:
        print("    Se encontraron errores de formato (muestra de hasta 2000 lineas):")
        for e in errors[:10]:
            print(f"      - {e}")
        raise RuntimeError(
            f"El BED tiene {len(errors)} problemas de formato detectados en el muestreo. "
            "Corrige el CSV de entrada o la logica de limpieza antes de continuar."
        )

    print("    Formato del BED verificado correctamente (chr*, enteros, 4 columnas).")


# --------------------------------------------------------------------------
# Paso 6: ejecucion de deepTools computeMatrix
# --------------------------------------------------------------------------
def run_compute_matrix(bed_path: Path, bigwig_path: Path, out_matrix: Path) -> None:
    print("[6/6] Ejecutando deepTools computeMatrix (reference-point, center)")

    if not bigwig_path.exists():
        raise FileNotFoundError(f"No se encontro el bigWig de ATAC-seq: {bigwig_path}")

    out_matrix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "computeMatrix", "reference-point",
        "-S", str(bigwig_path),
        "-R", str(bed_path),
        "--referencePoint", "center",
        "-b", "500",
        "-a", "500",
        "--skipZeros",
        "--numberOfProcessors", "4",
        "-o", str(out_matrix),
    ]

    print(f"    Comando: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "No se encontro el ejecutable 'computeMatrix'. Verifica que el ambiente "
            "'chromacrispr-phase1' este activo y que deepTools este instalado "
            "(conda activate chromacrispr-phase1)."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"computeMatrix fallo con codigo de salida {exc.returncode}. "
            "Revisa el mensaje de deepTools arriba para el detalle del error."
        ) from exc

    if not out_matrix.exists() or out_matrix.stat().st_size == 0:
        raise RuntimeError(
            f"computeMatrix termino sin error pero no genero un archivo valido en {out_matrix}"
        )

    print(f"    Matriz generada correctamente: {out_matrix}")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> None:
    print("=" * 70)
    print("ChromaCRISPR Phase 1 - F2: Preparacion ATAC-seq (deepTools)")
    print("=" * 70)

    try:
        df = load_and_filter(INPUT_CSV)
        df = clean_and_validate(df)
        bed = build_bed(df)
        write_bed(bed, OUTPUT_BED)
        validate_bed_file(OUTPUT_BED)
        run_compute_matrix(OUTPUT_BED, ATAC_BIGWIG, OUTPUT_MATRIX)
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print("\nProceso completado exitosamente.")
    print(f"  BED final:    {OUTPUT_BED}")
    print(f"  Matriz ATAC:  {OUTPUT_MATRIX}")


if __name__ == "__main__":
    main()
