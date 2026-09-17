#!/usr/bin/env python3

import pandas as pd
from pyfaidx import Fasta


def reverse_complement(seq):
    comp = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )
    return seq.translate(comp)[::-1]


def find_with_1_mismatch(guide, sequence):
    hits = []

    for i in range(len(sequence) - len(guide) + 1):
        sub = sequence[i:i+len(guide)]

        mismatches = sum(
            a != b
            for a, b in zip(guide, sub)
        )

        if mismatches <= 1:
            hits.append(
                (i, mismatches)
            )

    return hits


# inputs

unresolved_file = "logs/sanson_tss_resolution_complete.tsv"

tss_file = "data/reference/hg38/gencode_v49_TSS_nochr.tsv"

fasta_file = (
    "data/reference/hg38/fasta/"
    "Homo_sapiens.GRCh38.primary_assembly.fa"
)


df = pd.read_csv(
    unresolved_file,
    sep="\t"
)


df = df[
    df["status"]=="unresolved"
]


print("Unresolved actuales:", len(df))


tss = pd.read_csv(
    tss_file,
    sep="\t"
)


tss = (
    tss
    .drop_duplicates(
        subset=["gene_name"],
        keep=False
    )
)


df = df.merge(
    tss[
        [
            "gene_name",
            "chromosome",
            "tss"
        ]
    ],
    left_on="gene_symbol",
    right_on="gene_name",
    how="left"
)


fasta = Fasta(fasta_file)


WINDOW = 2000


resolved = []
failed = []


for _, row in df.iterrows():

    gene = row["gene_symbol"]
    guide = row["guide_sequence"]


    if pd.isna(row["tss"]):
        failed.append(
            {
                "gene_symbol":gene,
                "guide_sequence":guide,
                "reason":"no_gencode_tss"
            }
        )
        continue


    chrom = str(row["chromosome"])

    start = max(
        0,
        int(row["tss"]) - WINDOW
    )

    end = int(row["tss"]) + WINDOW


    seq = str(
        fasta[chrom][start:end]
    ).upper()


    searches = [
        (guide, "+"),
        (reverse_complement(guide), "-")
    ]


    found = False


    for target, strand in searches:

        hits = find_with_1_mismatch(
            guide,
            seq
        )


        if hits:

            pos, mismatches = hits[0]


            guide_pos = start + pos

            offset = (
                guide_pos -
                int(row["tss"])
            )


            resolved.append(
                {
                    "gene_symbol":gene,
                    "guide_sequence":guide,
                    "strand":strand,
                    "guide_position":guide_pos,
                    "offset":offset,
                    "mismatches":mismatches
                }
            )

            found=True
            break


    if not found:
        failed.append(
            {
                "gene_symbol":gene,
                "guide_sequence":guide,
                "reason":"no_match_1mm"
            }
        )


print("\n===== RESULT =====")
print("Recovered:", len(resolved))
print("Failed:", len(failed))


pd.DataFrame(resolved).to_csv(
    "logs/sanson_mismatch_resolved.tsv",
    sep="\t",
    index=False
)


pd.DataFrame(failed).to_csv(
    "logs/sanson_mismatch_failed.tsv",
    sep="\t",
    index=False
)
