#!/usr/bin/env python3

import pandas as pd
import argparse
from pyfaidx import Fasta


def reverse_complement(seq):
    complement = str.maketrans(
        "ACGT",
        "TGCA"
    )
    return seq.translate(complement)[::-1]


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--sgrna", required=True)
    parser.add_argument("--tss", required=True)
    parser.add_argument("--fasta", required=True)

    parser.add_argument(
        "--n",
        type=int,
        default=100
    )

    parser.add_argument(
        "--report",
        type=str,
        default=None,
        help="Ruta de salida TSV con clasificación de fallos de mapeo"
    )

    args = parser.parse_args()


    # Load sgRNA
    sgrna = pd.read_csv(
        args.sgrna,
        low_memory=False
    )

    sgrna = sgrna[
        sgrna["dataset"] == "Sanson2018"
    ].copy()


    # Test subset
    if args.n > 0:
        sgrna = sgrna.head(args.n)


    # Load TSS
    tss = pd.read_csv(
        args.tss,
        sep="\t"
    )

    # --- Detectar ambigüedad ANTES de deduplicar ---
    # Genes con más de un TSS anotado en GENCODE se resolvían en silencio
    # tomando el primero. Ahora se marcan explícitamente.
    tss_counts = tss.groupby("gene_name").size()
    ambiguous_genes = set(tss_counts[tss_counts > 1].index)

    tss_unique = (
        tss
        .drop_duplicates(
            subset=["gene_name"]
        )
    )


    merged = sgrna.drop(
        columns=["chromosome", "strand"],
        errors="ignore"
    ).merge(
        tss_unique[
            [
            "gene_name",
            "chromosome",
            "strand",
            "tss"
            ]
        ],
        left_on="gene_symbol",
        right_on="gene_name",
        how="left",
        indicator=True   # para distinguir "gene no encontrado" de "gene encontrado, tss nulo"
    )


    fasta = Fasta(args.fasta)


    mapped = 0
    results = []
    failed_records = []

    for _, row in merged.iterrows():

        gene = row["gene_symbol"]
        guide = row["guide_sequence"]

        # Caso 1: el gene_symbol no existe en absoluto en la tabla TSS
        if row["_merge"] == "left_only":
            failed_records.append({
                "guide_sequence": guide,
                "gene_symbol": gene,
                "failure_type": "missing_gene",
                "detail": "gene_symbol no encontrado en gencode_v49_TSS_nochr.tsv"
            })
            continue

        # Caso 2: el gen existía en múltiples TSS y se resolvió arbitrariamente
        if gene in ambiguous_genes:
            failed_records.append({
                "guide_sequence": guide,
                "gene_symbol": gene,
                "failure_type": "ambiguous_gene",
                "detail": f"{tss_counts[gene]} TSS candidatos en GENCODE para este gen"
            })
            continue

        # Caso 3: el gen matcheó pero el campo tss es nulo en la fuente
        if pd.isna(row["tss"]):
            failed_records.append({
                "guide_sequence": guide,
                "gene_symbol": gene,
                "failure_type": "missing_tss",
                "detail": "gene encontrado pero campo tss nulo en GENCODE"
            })
            continue

        chrom = str(row["chromosome"])
        tss_pos = int(row["tss"])
        strand = row["strand"]

        # window
        if strand == "+":
            start = max(0, tss_pos - 1)
            end = tss_pos + 300
        else:
            start = max(0, tss_pos - 300)
            end = tss_pos

        try:
            seq = str(
                fasta[chrom][start:end]
            ).upper()

        except Exception as e:
            # Caso 4: el cromosoma no existe en el FASTA o la ventana es inválida
            failed_records.append({
                "guide_sequence": guide,
                "gene_symbol": gene,
                "failure_type": "sequence_absent",
                "detail": f"error extrayendo ventana FASTA ({chrom}:{start}-{end}): {e}"
            })
            continue

        rc = reverse_complement(guide)

        found = None
        hit_strand = None

        if guide in seq:
            found = guide
            hit_strand = "+"
        elif rc in seq:
            found = rc
            hit_strand = "-"

        if found:
            mapped += 1
            results.append(
                {
                "guide_sequence": guide,
                "gene_symbol": gene,
                "chromosome": chrom,
                "coordinate_window_start": start,
                "coordinate_window_end": end,
                "tss": tss_pos,
                "strand": hit_strand
                }
            )
        else:
            # Caso 5: TSS válido, ventana extraída, pero la guía no aparece
            # en ninguna hebra dentro de la ventana de diseño
            failed_records.append({
                "guide_sequence": guide,
                "gene_symbol": gene,
                "failure_type": "sequence_absent",
                "detail": f"sin match en ventana {chrom}:{start}-{end} (ninguna hebra)"
            })


    print("Test guides:")
    print(len(merged))

    print("\nMapped:")
    print(mapped)

    print(
        "\nMapping rate:",
        mapped / len(merged) * 100,
        "%"
    )


    if results:
        print("\nExamples:")
        print(
            pd.DataFrame(results)
            .head(10)
            .to_string(index=False)
        )

        if results:
            mapped_df = pd.DataFrame(results)

            mapped_df.to_csv(
                "logs/sanson_window_mapped.tsv",
                sep="\t",
                index=False
            )

            print("\nReporte de mapeados escrito:")
            print("logs/sanson_window_mapped.tsv")

    if args.report:
        failed_df = pd.DataFrame(failed_records)
        failed_df.to_csv(
            args.report,
            sep="\t",
            index=False
        )

        print("\nReporte de fallos escrito:")
        print(args.report)

        if len(failed_df) > 0:
            print("\nTipos de fallo:")
            print(failed_df["failure_type"].value_counts())

if __name__ == "__main__":
    main()