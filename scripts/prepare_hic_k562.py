#!/usr/bin/env python3
"""
prepare_hic_k562.py

ChromaCRISPR Phase 1 - Semana 3

Preparación inicial Hi-C K562:
.hic -> .mcool -> extracción cool -> KR balancing

Objetivo:
SO1 Dataset curation and harmonization
Module 2 Epigenome Data Acquisition
"""

from pathlib import Path
import subprocess
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_HIC = (
    PROJECT_ROOT /
    "data/raw/hic_k562/4DNFITUOMFUQ.hic"
)

OUTPUT_DIR = (
    PROJECT_ROOT /
    "data/processed/hic_k562"
)


def run(cmd):

    print("\nRunning:")
    print(" ".join(cmd))

    result = subprocess.run(cmd)

    if result.returncode != 0:
        sys.exit(
            f"ERROR executing command: {' '.join(cmd)}"
        )


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not INPUT_HIC.exists():
        sys.exit(
            f"Missing Hi-C file:\n{INPUT_HIC}"
        )


    mcool_file = (
        OUTPUT_DIR /
        "4DNFITUOMFUQ.mcool"
    )

    cool_file = (
        OUTPUT_DIR /
        "4DNFITUOMFUQ.cool"
    )


    # Step 1
    # Convert hic to mcool

    if not mcool_file.exists():

        run([
            "hic2cool",
            "convert",
            str(INPUT_HIC),
            str(mcool_file)
        ])

    else:

        print(
            "Existing mcool detected."
        )


    # Step 2
    # Extract default cool resolution

    if not cool_file.exists():

        run([
            "cooler",
            "zoomify",
            str(mcool_file),
            "-o",
            str(cool_file)
        ])

    else:

        print(
            "Existing cool detected."
        )


    # Step 3
    # KR balancing

    run([
        "cooler",
        "balance",
        "--mad-max",
        "5",
        "-p",
        "4",
        str(cool_file)
    ])


    print(
        "\nHi-C preparation completed."
    )

    print(
        f"Output directory:\n{OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()
