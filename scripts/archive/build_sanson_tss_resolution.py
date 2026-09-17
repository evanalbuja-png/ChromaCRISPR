#!/usr/bin/env python3

import pandas as pd


# -----------------------------
# Inputs
# -----------------------------

gencode = pd.read_csv(
    "logs/window2kb_resolved_positions.tsv",
    sep="\t"
)

fantom = pd.read_csv(
    "logs/fantom5_resolved.tsv",
    sep="\t"
)


# -----------------------------
# GENCODE resolved
# -----------------------------

gencode["coordinate_source"] = "GENCODE"

gencode = gencode.rename(
    columns={
        "guide_position": "guide_position",
        "gencode_tss": "tss_position"
    }
)

gencode["crispri_region"] = (
    (gencode["offset"] >= -50) &
    (gencode["offset"] <= 300)
)

gencode["status"] = gencode["crispri_region"].map(
    {
        True: "promoter",
        False: "distal"
    }
)


gencode_final = gencode[
    [
        "gene_symbol",
        "guide_sequence",
        "coordinate_source",
        "guide_position",
        "tss_position",
        "offset",
        "crispri_region",
        "status"
    ]
]


# -----------------------------
# FANTOM5 resolved
# -----------------------------

fantom["coordinate_source"] = "FANTOM5"

fantom = fantom.rename(
    columns={
        "fantom_tss": "tss_position"
    }
)

fantom["crispri_region"] = (
    (fantom["offset"] >= -50) &
    (fantom["offset"] <= 300)
)

fantom["status"] = fantom["crispri_region"].map(
    {
        True: "promoter",
        False: "distal"
    }
)


fantom_final = fantom[
    [
        "gene_symbol",
        "guide_sequence",
        "coordinate_source",
        "guide_position",
        "tss_position",
        "offset",
        "crispri_region",
        "status"
    ]
]


# -----------------------------
# Merge
# -----------------------------

final = pd.concat(
    [
        gencode_final,
        fantom_final
    ],
    ignore_index=True
)


# eliminar duplicados exactos
final = final.drop_duplicates()


final.to_csv(
    "logs/sanson_tss_resolution_final.tsv",
    sep="\t",
    index=False
)


print("===== Final TSS Resolution =====")
print("Total resolved:", len(final))

print("\nPor fuente:")
print(final["coordinate_source"].value_counts())

print("\nPor categoría:")
print(final["status"].value_counts())

print("\nCRISPRi promoter:")
print(final["crispri_region"].sum())
