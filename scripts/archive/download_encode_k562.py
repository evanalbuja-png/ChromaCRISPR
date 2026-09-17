#!/usr/bin/env python3
"""
ChromaCRISPR Phase 1
Module 2 - ENCODE K562 acquisition

Downloads ENCODE bigWig signal tracks using ENCODE accessions.
"""

from pathlib import Path
import requests
import logging

# ---------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------

OUTPUT_DIR = Path("data/raw/encode_k562")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LOGFILE = OUTPUT_DIR / "download_log.txt"

logging.basicConfig(
    filename=LOGFILE,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)

# ---------------------------------------------------------------------
# ENCODE experiments from proposal
# ---------------------------------------------------------------------

ACCESSIONS = {
    "ATAC": "ENCSR868FGK",
    "H3K27ac": "ENCSR000AKP",
    "H3K4me3": "ENCSR000AKQ",
    "H3K27me3": "ENCSR000AKS",
    "H3K9me3": "ENCSR000DWD",
}

HEADERS = {
    "accept": "application/json"
}

# ---------------------------------------------------------------------
# ENCODE API
# ---------------------------------------------------------------------

def get_experiment_metadata(accession):
    """
    Query ENCODE API for an experiment accession.

    Returns
    -------
    dict
        JSON metadata describing the experiment.
    """

    url = (
        f"https://www.encodeproject.org/experiments/"
        f"{accession}/?format=json"
    )

    print(f"Querying {accession}...")

    response = requests.get(
        url,
        headers=HEADERS,
        timeout=60
    )

    response.raise_for_status()

    metadata = response.json()

    logging.info(f"Retrieved metadata for {accession}")

    return metadata

    # ---------------------------------------------------------------------
# Find the correct bigWig signal file
# ---------------------------------------------------------------------

def find_bigwig(metadata):
    """
    Locate the preferred bigWig file within an ENCODE experiment.

    Preference:
        - file_format == "bigWig"
        - status == "released"

    Returns
    -------
    dict
        Metadata for the selected file.
    """

    files = metadata.get("files", [])

    for f in files:

        if (
            f.get("file_format") == "bigWig"
            and f.get("status") == "released"
        ):
            return f

    raise RuntimeError(
        f"No released bigWig found for "
        f"{metadata.get('accession')}"
    )

    # ---------------------------------------------------------------------
# Download a file from ENCODE
# ---------------------------------------------------------------------

def download_file(file_info):
    """
    Download a released ENCODE file.

    Parameters
    ----------
    file_info : dict
        Dictionary returned by find_bigwig().
    """

    href = file_info["href"]

    download_url = f"https://www.encodeproject.org{href}"

    filename = href.split("/")[-1]

    outfile = OUTPUT_DIR / filename

    if outfile.exists():

        print(f"✓ {filename} already exists.")

        logging.info(f"Skipped {filename}")

        return outfile

    print(f"Downloading {filename}...")

    response = requests.get(
        download_url,
        stream=True,
        timeout=120
    )

    response.raise_for_status()

    with open(outfile, "wb") as f:

        for chunk in response.iter_content(chunk_size=1024 * 1024):

            if chunk:

                f.write(chunk)

    size_mb = outfile.stat().st_size / (1024 * 1024)

    print(f"Finished ({size_mb:.2f} MB)")

    logging.info(
        f"Downloaded {filename} "
        f"({size_mb:.2f} MB)"
    )

    return outfile

    # ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print("=" * 60)
    print("ChromaCRISPR Phase 1")
    print("ENCODE K562 Acquisition")
    print("=" * 60)

    logging.info("Starting ENCODE download")

    for assay, accession in ACCESSIONS.items():

        print(f"\nProcessing {assay} ({accession})")

        try:

            metadata = get_experiment_metadata(accession)

            bigwig = find_bigwig(metadata)

            outfile = download_file(bigwig)

            print(f"Saved to: {outfile}")

        except Exception as e:

            print(f"ERROR: {assay} -> {e}")

            logging.exception(
                f"Failed downloading {assay} ({accession})"
            )

    print("\nFinished.")
    logging.info("Download finished")


if __name__ == "__main__":
    main()