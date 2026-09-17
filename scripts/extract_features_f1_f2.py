#!/usr/bin/env python3

"""
ChromaCRISPR Phase 1
Week 5
Feature extraction:
F1 - sequence features
F2 - ATAC accessibility

Input:
data/interim/sgRNA_unified.csv

Output:
data/interim/features_f1_f2.h5
"""


import os
import subprocess
import pandas as pd
import numpy as np
import h5py


INPUT = "data/interim/sgRNA_unified.csv"

OUTPUT = "data/interim/features_f1_f2.h5"

BED_FILE = "data/interim/features/sgRNA_hg38.bed"

ATAC_BIGWIG = (
    "data/raw/encode_k562/"
    "ATAC-seq_K562.bigWig"
)


WINDOWS = [
    100,
    200,
    500
]


def gc_content(seq):

    seq = seq.upper()

    return (
        seq.count("G") +
        seq.count("C")
    ) / len(seq)



def sequence_features(df):

    print("Extracting F1 sequence features")

    df["gc_content"] = (
        df["guide_sequence"]
        .apply(gc_content)
    )

    df["guide_length_calc"] = (
        df["guide_sequence"]
        .str.len()
    )

    df["pam_GG"] = (
        df["guide_sequence"]
        .str.endswith("GG")
        .astype(int)
    )

    df["pam_GC"] = (
        df["guide_sequence"]
        .str[-2:]
        .isin(["GC","CG"])
        .astype(int)
    )

    # NUEVAS FEATURES TSS

    df["tss_offset"] = df["offset"]

    df["crispri_region"] = (
        df["crispri_region"]
        .astype(int)
    )

    df["tss_source_encoded"] = (
        df["coordinate_source"]
        .map({
            "GENCODE":0,
            "FANTOM5":1,
            "unresolved_TSS_review":2
        })
    )

    return df



def create_hg38_bed(df):

    print("Creating hg38 BED")

    hg38 = df[
        df.genome_build=="hg38"
    ].copy()


    hg38 = hg38.dropna(
        subset=[
            "chromosome",
            "coordinate"
        ]
    )


    bed = pd.DataFrame({

        "chrom":
        hg38.chromosome,

        "start":
        hg38.coordinate.astype(int)-1,

        "end":
        hg38.coordinate.astype(int),

        "name":
        hg38.guide_sequence

    })


    bed.to_csv(
        BED_FILE,
        sep="\t",
        index=False,
        header=False
    )


    return hg38



def run_deeptools():

    print(
        "Running deepTools computeMatrix"
    )

    output = (
        "data/interim/features/"
        "atac_matrix.gz"
    )


    cmd = [
        "computeMatrix",
        "reference-point",
        "-S",
        ATAC_BIGWIG,
        "-R",
        BED_FILE,
        "--referencePoint",
        "center",
        "-b",
        "500",
        "-a",
        "500",
        "--skipZeros",
        "-o",
        output
    ]


    subprocess.run(
        cmd,
        check=True
    )


    return output



def main():

    print("Loading sgRNA dataset")

    df = pd.read_csv(
        INPUT,
        low_memory=False
    )


    print(
        f"Initial guides: {len(df)}"
    )


    df = sequence_features(df)


    hg38 = create_hg38_bed(df)


    atac_matrix = run_deeptools()


    print(
        "Saving feature matrix"
    )


    with h5py.File(
        OUTPUT,
        "w"
    ) as h5:

        h5.create_dataset(
            "gc_content",
            data=df.gc_content.values
        )


        h5.create_dataset(
            "guide_length",
            data=df.guide_length.values
        )

        h5.create_dataset(
            "tss_offset",
            data=df.tss_offset.values
        )

        h5.create_dataset(
            "crispri_region",
            data=df.crispri_region.values
        )

        h5.create_dataset(
            "tss_source_encoded",
            data=df.tss_source_encoded.values
        )

        h5.create_dataset(
            "hg38_guides",
            data=hg38.index.values
        )


    print(
        "DONE"
    )


if __name__=="__main__":

    main()
