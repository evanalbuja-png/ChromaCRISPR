#!/usr/bin/env python3
"""
prepare_hic_k562.py

Semana 3 - ChromaCRISPR Phase 1

Preparación inicial de datos Hi-C K562:
    .hic -> .cool -> cooler balance

Entrada:
    data/raw/hic_k562/4DNFITUOMFUQ.hic

Salida:
    data/processed/hic_k562/4DNFITUOMFUQ.cool
    data/processed/hic_k562/4DNFITUOMFUQ_balanced.cool

Requisitos:
    - hic2cool
    - cooler
"""

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INPUT = (
    PROJECT_ROOT /
    "data/raw/hic_k562/4DNFITUOMFUQ.hic"
)

DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT /
    "data/processed/hic_k562"
)


def run_command(command):
    """Execute shell command and stop on failure."""
    print("\nRunning:")
    print(" ".join(command))

    result = subprocess.run(command)

    if result.returncode != 0:
        sys.exit(
            f"\nERROR: Command failed with exit code {result.returncode}"
        )


def main():

    parser = argparse.ArgumentParser(
        description="Convert Hi-C .hic file to .cool and balance."
    )

    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
        help="Input .hic file"
    )

    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory"
    )

    args = parser.parse_args()


    hic_file = Path(args.input)
    output_dir = Path(args.output_dir)

    if not hic_file.exists():
        sys.exit(
            f"ERROR: Hi-C file not found:\n{hic_file}"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    cool_file = (
        output_dir /
        f"{hic_file.stem}.cool"
    )

    balanced_file = (
        output_dir /
        f"{hic_file.stem}_balanced.cool"
    )


    print("=" * 70)
    print("ChromaCRISPR Phase 1 - Hi-C preparation")
    print("=" * 70)

    print(f"Input Hi-C:")
    print(hic_file)

    print(f"\nOutput directory:")
    print(output_dir)


    # ------------------------------------------------------
    # Step 1: Convert .hic -> .cool
    # ------------------------------------------------------

    if not cool_file.exists():

        run_command([
            "hic2cool",
            "convert",
            str(hic_file),
            str(cool_file)
        ])

    else:
        print(
            f"\nExisting COOL file found:"
            f"\n{cool_file}"
        )


    # ------------------------------------------------------
    # Step 2: KR balancing
    # ------------------------------------------------------

    if not balanced_file.exists():

        run_command([
            "cooler",
            "balance",
            "--mad-max",
            "5",
            "-p",
            "4",
            str(cool_file)
        ])

        # cooler balance modifies the file in place.
        # Create a copy with explicit balanced name.

        run_command([
            "cp",
            str(cool_file),
            str(balanced_file)
        ])

    else:
        print(
            f"\nExisting balanced COOL found:"
            f"\n{balanced_file}"
        )


    print("\n" + "=" * 70)
    print("Hi-C preparation completed")
    print("=" * 70)

    print("\nGenerated files:")

    for f in [
        cool_file,
        balanced_file
    ]:
        if f.exists():
            print(
                f"{f} "
                f"({f.stat().st_size / 1024**3:.2f} GB)"
            )


if __name__ == "__main__":
    main()
