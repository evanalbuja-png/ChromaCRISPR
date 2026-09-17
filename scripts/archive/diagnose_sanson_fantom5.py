#!/usr/bin/env python3

import pandas as pd
from pyfaidx import Fasta


def reverse_complement(seq):
    comp = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )
    return seq.translate(comp)[::-1]


# ==========================
# Inputs
# ==========================

unresolved_file = "logs/window2kb_unresolved.tsv"

fantom_file = (
    "data/reference/fantom5/fantom5_peaks.tsv"
)

fasta_file = (
    "data/reference/hg38/fasta/"
    "Homo_sapiens.GRCh38.primary_assembly.fa"
)


# ==========================
# Load data
# ==========================

unresolved = pd.read_csv(
    unresolved_file,
    sep="\t"
)

fantom = pd.read_csv(
    fantom_file,
    sep="\t"
)


fasta = Fasta(fasta_file)


print("Unresolved:", len(unresolved))
print("FANTOM peaks:", len(fantom))


# ==========================
# Parameters
# ==========================

WINDOW = 2000


resolved = []
failed = []


# ==========================
# Search
# ==========================

for _, row in unresolved.iterrows():

    gene = row["gene_symbol"]
    guide = row["guide_sequence"]

    rc = reverse_complement(guide)


    gene_peaks = fantom[
        fantom["gene_symbol"] == gene
    ]


    if len(gene_peaks) == 0:

        failed.append({
            "gene_symbol": gene,
            "guide_sequence": guide,
            "reason": "no_fantom_peak"
        })

        continue


    found = False


    for _, peak in gene_peaks.iterrows():

        chrom = str(peak["chromosome"]).replace("chr", "")
        fantom_tss = int(
            peak["fantom_tss"]
        )


        start = max(
            0,
            fantom_tss - WINDOW
        )

        end = fantom_tss + WINDOW


        try:

            seq = str(
                fasta[chrom][start:end]
            ).upper()

        except Exception:

            continue


        pos_fwd = seq.find(guide)
        pos_rev = seq.find(rc)


        if pos_fwd != -1:

            guide_pos = start + pos_fwd

            if peak["strand"] == "+":
                offset = guide_pos - fantom_tss
            else:
                offset = fantom_tss - guide_pos

            resolved.append({

                "gene_symbol": gene,
                "guide_sequence": guide,
                "strand_found": "+",
                "guide_position": guide_pos,
                "fantom_tss": fantom_tss,
                "offset": offset,
                "peak_score": peak["score"]

            })

            found = True
            break


        elif pos_rev != -1:

            guide_pos = start + pos_rev

            if peak["strand"] == "+":
                offset = guide_pos - fantom_tss
            else:
                offset = fantom_tss - guide_pos

            resolved.append({

                "gene_symbol": gene,
                "guide_sequence": guide,
                "strand_found": "-",
                "guide_position": guide_pos,
                "fantom_tss": fantom_tss,
                "offset": offset,
                "peak_score": peak["score"]

            })

            found = True
            break


    if not found:

        failed.append({
            "gene_symbol": gene,
            "guide_sequence": guide,
            "reason": "no_match_fantom5",
            "fantom_peaks_tested": len(gene_peaks)
        })

        if len(failed) <= 10:
            print(
                "NO MATCH:",
                gene,
                "peaks:",
                len(gene_peaks)
            )


# ==========================
# Output
# ==========================

print("\n===== Summary =====")
print("Resolved:", len(resolved))
print("Failed:", len(failed))


pd.DataFrame(resolved).to_csv(
    "logs/fantom5_resolved.tsv",
    sep="\t",
    index=False
)


pd.DataFrame(failed).to_csv(
    "logs/fantom5_unresolved.tsv",
    sep="\t",
    index=False
)