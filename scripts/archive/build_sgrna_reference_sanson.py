#!/usr/bin/env python3
"""
build_sgrna_reference_FINAL.py
Reconstruye la referencia genómica de sgRNA desde cero, en una sola pasada,
con todas las reglas validadas durante Semana 5.
Uso: python build_sgrna_reference_FINAL.py --dataset Sanson2018
     python build_sgrna_reference_FINAL.py --dataset Replogle2022
"""

import pandas as pd
import argparse
from pyfaidx import Fasta

# ============================================================
# REGLAS FIJAS (validadas hoy — no se vuelven a discutir)
# ============================================================
WINDOW_UPSTREAM = 50    # -50pb, confirmado contra literatura CRISPRi
WINDOW_DOWNSTREAM = 300 # +300pb, confirmado contra literatura CRISPRi

# Familias de genes conocidas por fallar en anotación TSS
# (confirmado hoy: OR = expresión no capturada por FANTOM5;
#  ZNF/MUC/HLA = regiones repetitivas/polimórficas de mapeo ambiguo)
KNOWN_HARD_FAMILIES = ("OR", "ZNF", "MUC", "HLA", "HIST")

def strand_aware_offset(guide_pos, tss_pos, strand):
    """Offset corregido por strand desde el inicio — no se calcula nunca sin esto."""
    if strand == "+":
        return guide_pos - tss_pos
    else:
        return tss_pos - guide_pos

def classify(offset):
    if offset is None:
        return "unresolved"
    return "promoter" if (-WINDOW_UPSTREAM <= offset <= WINDOW_DOWNSTREAM) else "distal"

def exact_match_in_window(guide_seq, chrom, center_pos, upstream, downstream, fasta):
    """Búsqueda exacta de substring, ambas hebras. NO usa Bowtie2 (evita
    ambigüedad de parámetros de mismatch en secuencias de 20nt)."""
    start = max(0, center_pos - upstream)
    end = center_pos + downstream
    try:
        fasta_chrom = chrom.replace("chr", "")
        window_seq = str(fasta[fasta_chrom][start:end]).upper()
    except Exception:
        return None, None, None

    rc = guide_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    pos = window_seq.find(guide_seq)
    if pos != -1:
        return start + pos, "+", "exact_window_fwd"
    pos = window_seq.find(rc)
    if pos != -1:
        return start + pos, "-", "exact_window_rev"
    return None, None, None

def resolve_guide(row, gencode_tss, fantom5_tss, fasta):
    gene = row["gene_symbol"]
    guide = row["guide_sequence"]

    if gene == "GATA3":
        print("DEBUG GATA3")
        print("guide:", guide)
        print("GENCODE:", gencode_tss.get(gene))
        print("FANTOM5:", fantom5_tss.get(gene))

    if gene == "GATA3":
        tss_pos, chrom, strand = fantom5_tss[gene]
        start = max(0, tss_pos - 50)
        end = tss_pos + 300
        print("FASTA chrom exists:", chrom in fasta)
        print("Window:", start, end)
        fasta_chrom = chrom.replace("chr", "")
        print("FASTA chrom:", fasta_chrom)
        print("FASTA exists:", fasta_chrom in fasta)
        print("FASTA length:", len(fasta[fasta_chrom]))
        print("Window seq:", str(fasta[fasta_chrom][start:end]).upper())
        print("Guide:", guide)

    # PASO 0: exclusión automática de familias conocidas como difíciles
    # (esto evita re-litigar OR7E24 / ZNF404 / MUC22 uno por uno otra vez)
    if any(gene.startswith(fam) for fam in KNOWN_HARD_FAMILIES):
        return {
            "coordinate_source": "excluded_hard_family",
            "status": "unresolved",
            "detail": f"Familia {gene[:3]} — conocida por anotación TSS no confiable"
        }

    # PASO 1: ventana estrecha contra TSS de GENCODE (rápido, cubre ~50%)
    if gene in gencode_tss:
        tss_pos, chrom, strand = gencode_tss[gene]
        pos, hit_strand, method = exact_match_in_window(
            guide, chrom, tss_pos, WINDOW_UPSTREAM, WINDOW_DOWNSTREAM, fasta
        )
        if pos is not None:
            offset = strand_aware_offset(pos, tss_pos, strand)
            return {
                "coordinate_source": "GENCODE_window",
                "chromosome": chrom, "coordinate": pos, "strand": strand,
                "offset": offset, "status": classify(offset)
            }

    # PASO 2: fallback a pico CAGE de FANTOM5 (cubre promotores alternativos)
    if gene in fantom5_tss:
        tss_pos, chrom, strand = fantom5_tss[gene]
        pos, hit_strand, method = exact_match_in_window(
            guide, chrom, tss_pos, WINDOW_UPSTREAM, WINDOW_DOWNSTREAM, fasta
        )
        if pos is not None:
            offset = strand_aware_offset(pos, tss_pos, strand)
            return {
                "coordinate_source": "FANTOM5_window",
                "chromosome": chrom, "coordinate": pos, "strand": strand,
                "offset": offset, "status": classify(offset)
            }

    # PASO 3: ventana ampliada (±2kb) como último recurso, MISMO gen únicamente
    # (nunca genoma completo — así se evita el multi-mapping ambiguo de Bowtie2)
    for tss_source_name, tss_source in [("GENCODE", gencode_tss), ("FANTOM5", fantom5_tss)]:
        if gene in tss_source:
            tss_pos, chrom, strand = tss_source[gene]
            pos, hit_strand, method = exact_match_in_window(
                guide, chrom, tss_pos, 2000, 2000, fasta
            )
            if pos is not None:
                offset = strand_aware_offset(pos, tss_pos, strand)
                return {
                    "coordinate_source": f"{tss_source_name}_wide_2kb",
                    "chromosome": chrom, "coordinate": pos, "strand": strand,
                    "offset": offset, "status": classify(offset)
                }

    # PASO 4: no se resolvió por ningún método legítimo — se documenta, no se inventa
    return {
        "coordinate_source": "unresolved_TSS_review",
        "status": "unresolved",
        "detail": f"{gene}: sin match en GENCODE/FANTOM5, ventana estrecha ni ±2kb"
    }

def validate_output(df, dataset_name):
    errors = []

    # 1. No debe haber duplicados guía+gen
    if df.duplicated(["guide_sequence", "gene_symbol"]).sum() > 0:
        errors.append("Duplicados guía+gen detectados")

    # 2. La mayoría de GENCODE_window debe ser 'promoter' (si no, algo estructural falló)
    window_only = df[df["coordinate_source"] == "GENCODE_window"]
    if len(window_only) > 0:
        promoter_rate = (window_only["status"] == "promoter").mean()
        if promoter_rate < 0.85:
            errors.append(f"ALERTA: solo {promoter_rate:.1%} de GENCODE_window es "
                          f"'promoter' — se esperaba >85%, revisar offset/strand")

    # 3. Tasa de resolución total razonable
    resolved_rate = (df["status"] != "unresolved").mean()
    print(f"Tasa de resolución {dataset_name}: {resolved_rate:.2%}")
    if resolved_rate < 0.85:
        errors.append(f"Tasa de resolución baja ({resolved_rate:.2%}) — revisar antes de continuar")

    if errors:
        print("\n🚩 VALIDACIÓN FALLÓ:")
        for e in errors:
            print(f"  - {e}")
        raise AssertionError("Corrige antes de continuar — no avances con datos no validados")
    else:
        print("✅ Todas las validaciones pasaron")

    return df    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()

    if args.dataset != "Sanson2018":
        raise ValueError("Dataset no implementado todavía")

    print("=== Construyendo referencia Sanson2018 ===")

    # Cargar guía original
    guides = pd.read_csv(
        "data/interim/sanson_features_reference.tsv",
        sep="\t"
    )

    guides = guides[[
    "guide_sequence",
    "gene_symbol"
    ]]

    # Cargar GENCODE TSS
    gencode_df = pd.read_csv(
        "data/reference/hg38/gencode_v49_TSS.tsv",
        sep="\t"
    )

    gencode_tss = {
        row["gene_name"]: (
            row["tss"],
            row["chromosome"],
            row["strand"]
        )
        for _, row in gencode_df.iterrows()
    }

    print("GENCODE genes cargados:", len(gencode_tss))
    print("Ejemplo GENCODE:", list(gencode_tss.items())[:3])

    # Cargar FANTOM5
    fantom_df = pd.read_csv(
        "data/reference/fantom5/fantom5_lookup.tsv",
        sep="\t"
    )

    fantom5_tss = {
        row["gene_symbol"]: (
            row["fantom_tss"],
            row["chromosome"],
            row["strand"]
        )
        for _, row in fantom_df.iterrows()
    }

    print("FANTOM5 genes cargados:", len(fantom5_tss))
    print("Ejemplo FANTOM5:", list(fantom5_tss.items())[:3])

    # Cargar FASTA hg38
    fasta = Fasta(
        "data/reference/hg38/fasta/Homo_sapiens.GRCh38.primary_assembly.fa"
    )

    print(f"Guías a procesar: {len(guides)}")

    results = []

    for _, row in guides.iterrows():
        result = resolve_guide(
            row,
            gencode_tss,
            fantom5_tss,
            fasta
        )
        results.append(result)

    result_df = pd.DataFrame(results)

    final_df = pd.concat(
        [
            guides.reset_index(drop=True),
            result_df.reset_index(drop=True)
        ],
        axis=1
    )

    print(final_df["coordinate_source"].value_counts(dropna=False).head(20))
    print(final_df["status"].value_counts(dropna=False))
    print(final_df.head())

    validate_output(
        final_df,
        args.dataset
    )

    output = "data/processed/sanson2018_reference_FINAL.csv"

    final_df.to_csv(
        output,
        index=False
    )

    print(f"Archivo generado: {output}")