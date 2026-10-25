from importlib.metadata import version

from pathlib import Path
from datetime import datetime, timezone
import csv
import json
import sys

import scanpy as sc

# ==========================================
# MI-Agata-MultiOmicsFlow
# Month 3 - Single-cell RNA-seq
# Step 01: Dataset Assessment
# Dataset: GEO GSE176078 (Wu et al., 2021)
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
REPORT_DIR = BASE_DIR / "reports"

DATASET_ID = "GSE176078"
DATASET_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/"
    "query/acc.cgi?acc=GSE176078"
)

SUPPORTED = {".h5ad", ".h5", ".mtx", ".gz", ".csv", ".tsv"}

def main():
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n=== SINGLE-CELL DATASET ASSESSMENT ===")
    print(f"Dataset: {DATASET_ID}")
    print(f"Scanpy version: {version("scanpy")}")
    print(f"Python version: {sys.version.split()[0]}")
    print(f"Data directory: {DATA_DIR}")

    files = sorted(
        (p for p in DATA_DIR.rglob("*") if p.is_file()),
        key=lambda p: str(p)
    ) if DATA_DIR.exists() else []

    rows = []

    for file in files:
        rows.append({
            "file": file.relative_to(DATA_DIR).as_posix(),
            "size_mb": round(file.stat().st_size / 1048576, 3),
            "supported_extension": file.suffix.lower() in SUPPORTED
        })

    print(f"\nFiles found: {len(files)}")

    for row in rows:
        print(f" - {row['file']} ({row['size_mb']} MB)")

    status = "FILES_FOUND" if files else "AWAITING_DATA"

    report = {
        "project": "MI-Agata-MultiOmicsFlow",
        "module": "Month 3 - Single-cell RNA-seq",
        "dataset_id": DATASET_ID,
        "dataset_url": DATASET_URL,
        "assessment_date_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "scanpy_version": version("scanpy"),
        "python_version": sys.version.split()[0],
        "files_found": len(files),
        "total_size_mb": round(
            sum(p.stat().st_size for p in files) / 1048576, 3
        ),
        "status": status,
        "note": (
            "Inventory only. Expression matrix and metadata "
            "have not yet been validated."
        )
    }

    report_path = REPORT_DIR / "dataset_assessment.json"
    inventory_path = REPORT_DIR / "dataset_file_inventory.csv"

    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    with inventory_path.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["file", "size_mb", "supported_extension"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nStatus: {status}")
    print(f"Report: {report_path}")
    print(f"Inventory: {inventory_path}")
    print("\n=== ASSESSMENT FINISHED ===")


if __name__ == "__main__":
    main()

