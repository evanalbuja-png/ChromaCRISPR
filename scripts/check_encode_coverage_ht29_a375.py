import requests

BASE_URL = "https://www.encodeproject.org/search/"

CELL_LINES = ["HT-29", "A375"]

MARKS = [
    "H3K27ac",
    "H3K4me3",
    "H3K4me1",
    "H3K27me3",
    "H3K9me3",
]

for cell in CELL_LINES:
    print("\n" + "=" * 70)
    print(f"CELL LINE: {cell}")
    print("=" * 70)

    params = {
        "type": "Experiment",
        "biosample_ontology.term_name": cell,
        "format": "json",
        "limit": 100,
    }

    response = requests.get(
        BASE_URL,
        params=params,
        headers={"Accept": "application/json"},
        timeout=30,
    )

    print("HTTP:", response.status_code)

    data = response.json()

    print("Total experiments:", data.get("total"))

    for experiment in data.get("@graph", []):
        print(
            experiment.get("accession"),
            "|",
            experiment.get("status"),
            "|",
            experiment.get("assay_title"),
        )

    print("\nHistone marks checked:")

    for mark in MARKS:
        print(f"  - {mark}")
