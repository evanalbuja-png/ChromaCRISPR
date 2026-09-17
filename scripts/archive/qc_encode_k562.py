#!/usr/bin/env python3
"""
qc_encode_k562.py

Quality control for ENCODE K562 BigWig files.

Checks:
    - File integrity (can pyBigWig open it?)
    - Header information
    - Chromosome names and lengths
    - Coverage statistics
    - Approximate min/max signal
    - Global mean signal

Output:
    data/interim/qc_encode_report.md
"""

from pathlib import Path
import math
import statistics

import pyBigWig


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

BIGWIG_DIR = PROJECT_ROOT / "data/raw/encode_k562"

REPORT_DIR = PROJECT_ROOT / "data/interim"

REPORT_FILE = REPORT_DIR / "qc_encode_report.md"


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def estimate_signal_range(bw, chrom_sizes, sample_window=10000):
    """
    Estimate min/max signal by sampling each chromosome.

    Reading an entire 1.6 GB bigWig would be unnecessarily expensive.
    This function samples one window near the beginning of each chromosome.
    """

    observed = []

    for chrom, length in chrom_sizes.items():

        end = min(sample_window, length)

        try:
            values = bw.values(chrom, 0, end)

            values = [
                x for x in values
                if x is not None and not math.isnan(x)
            ]

            observed.extend(values)

        except RuntimeError:
            continue

    if not observed:
        return None, None

    return min(observed), max(observed)


def compute_global_mean(header):
    """
    Compute global mean signal using header statistics.
    """

    covered = header.get("nBasesCovered", 0)

    if covered == 0:
        return None

    return header["sumData"] / covered


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    bigwigs = sorted(BIGWIG_DIR.glob("*.bigWig"))

    if not bigwigs:
        raise FileNotFoundError("No BigWig files found.")

    report = []

    report.append("# ENCODE K562 Quality Control Report\n")
    report.append("\n")
    report.append("Generated automatically by `qc_encode_k562.py`\n")
    report.append("\n")
    report.append("---\n")

    for bw_path in bigwigs:

        report.append(f"\n# {bw_path.name}\n")

        try:

            bw = pyBigWig.open(str(bw_path))

        except Exception as e:

            report.append(f"❌ Could not open file.\n")
            report.append(f"\nError:\n\n```\n{e}\n```\n")
            continue

        header = bw.header()

        chroms = bw.chroms()

        mean_signal = compute_global_mean(header)

        min_signal, max_signal = estimate_signal_range(
            bw,
            chroms
        )

        report.append("## Status\n")
        report.append("✅ File opened successfully.\n\n")

        report.append("## Header\n")

        for key, value in header.items():
            report.append(f"- **{key}**: {value}\n")

        report.append("\n")

        report.append("## Chromosomes\n")
        report.append(f"- Total chromosomes: {len(chroms)}\n")

        report.append("\n")

        report.append("| Chromosome | Length |\n")
        report.append("|------------|-------:|\n")

        for chrom, length in chroms.items():

            report.append(
                f"| {chrom} | {length:,} |\n"
            )

        report.append("\n")

        report.append("## Signal summary\n")

        report.append(
            f"- Approximate minimum signal: "
            f"{min_signal:.4f}\n"
            if min_signal is not None else
            "- Approximate minimum signal: N/A\n"
        )

        report.append(
            f"- Approximate maximum signal: "
            f"{max_signal:.4f}\n"
            if max_signal is not None else
            "- Approximate maximum signal: N/A\n"
        )

        report.append(
            f"- Global mean signal: "
            f"{mean_signal:.4f}\n"
            if mean_signal is not None else
            "- Global mean signal: N/A\n"
        )

        report.append(
            f"- Bases covered: "
            f"{header['nBasesCovered']:,}\n"
        )

        report.append("\n---\n")

        bw.close()

    REPORT_FILE.write_text(
        "".join(report),
        encoding="utf-8"
    )

    print("=" * 60)
    print("QC finished successfully.")
    print(f"Report written to:\n{REPORT_FILE}")
    print("=" * 60)


if __name__ == "__main__":
    main()