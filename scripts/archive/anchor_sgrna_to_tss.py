#!/usr/bin/env python3

import pandas as pd
import argparse


WINDOWS = {
    "Sanson2018": {
        "upstream": 0,
        "downstream": 300
    },
    "Replogle2022": {
        "upstream": 25,
        "downstream": 500
    }
}


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--dataset", required=True)
    parser.add_argument("--sgrna", required=True)
    parser.add_argument("--tss", required=True)

    args = parser.parse_args()

    sgrna = pd.read_csv(
        args.sgrna,
        low_memory=False
    )

    tss = pd.read_csv(
        args.tss,
        sep="\t"
    )

    print("Original sgRNAs:", len(sgrna))

    # dataset filter
    sgrna = sgrna[
        sgrna["dataset"] == args.dataset
    ].copy()

    print("Dataset rows:", len(sgrna))


    # merge gene symbol -> TSS
    # seleccionar un único TSS por gen
    tss_unique = (
        tss
        .sort_values("gene_id")
        .drop_duplicates(
            subset=["gene_name"],
            keep="first"
        )
    )


    merged = sgrna.merge(
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
        suffixes=("_old", "")
    )


    print("\nAfter TSS mapping:")
    print(
        merged["tss"].notna().sum(),
        "/",
        len(merged)
    )


    print("\nMissing TSS:")
    print(
        merged["tss"].isna().mean()*100,
        "%"
    )


    print("\nExample:")
    print(
        merged[
            [
            "guide_sequence",
            "gene_symbol",
            "chromosome",
            "strand",
            "tss"
            ]
        ].head(20)
    )


if __name__ == "__main__":
    main()
