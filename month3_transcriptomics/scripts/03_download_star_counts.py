
import hashlib
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data/raw/tcga_brca_selected_226.csv"
OUTPUT = ROOT / "data/raw/star_counts"

OUTPUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(METADATA, dtype=str)

required = {"file_id", "patient_id", "sample_type"}
assert required.issubset(df.columns)
assert len(df) == 226
assert df["file_id"].nunique() == 226

session = requests.Session()

def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def download_file(row):
    file_id = row["file_id"]
    patient = row["patient_id"]
    sample_type = row["sample_type"].replace(" ", "_")

    destination = OUTPUT / f"{patient}_{sample_type}_{file_id}.tsv"
    temporary = destination.with_suffix(".part")

    url = f"https://api.gdc.cancer.gov/data/{file_id}"

    expected_size = None
    if pd.notna(row.get("file_size_mb")):
        expected_size = float(row["file_size_mb"])

    for attempt in range(1, 4):
        try:
            if destination.exists():
                if expected_size is not None:
                    actual_mb = destination.stat().st_size / 1024**2
                    if abs(actual_mb - expected_size) < 0.02:
                        print(f"SKIP: {destination.name}", flush=True)
                        return "skipped"

            print(
                f"Downloading {patient} | {sample_type} "
                f"(attempt {attempt}/3)",
                flush=True
            )

            with session.get(
                url, stream=True, timeout=(30, 120)
            ) as response:
                response.raise_for_status()

                with open(temporary, "wb") as f:
                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):
                        if chunk:
                            f.write(chunk)

            if temporary.stat().st_size == 0:
                raise ValueError("Empty download")

            if expected_size is not None:
                actual_mb = temporary.stat().st_size / 1024**2
                if abs(actual_mb - expected_size) >= 0.02:
                    raise ValueError("File size mismatch")

            temporary.replace(destination)
            print(f"OK: {destination.name}", flush=True)
            return "downloaded"

        except (requests.RequestException, OSError, ValueError) as exc:
            print(f"ERROR: {exc}", flush=True)
            temporary.unlink(missing_ok=True)

            if attempt < 3:
                time.sleep(5 * attempt)

    return "failed"


def main():
    print("TCGA-BRCA STAR Counts download")
    print(f"Files: {len(df)}")

    results = []

    for i, (_, row) in enumerate(df.iterrows(), 1):
        print(f"\n[{i}/{len(df)}]", flush=True)

        status = download_file(row)
        results.append({
            "file_id": row["file_id"],
            "patient_id": row["patient_id"],
            "sample_type": row["sample_type"],
            "status": status
        })

    report = pd.DataFrame(results)
    report.to_csv(
        ROOT / "data/raw/star_counts_download_report.csv",
        index=False
    )

    print("\nDOWNLOAD SUMMARY")
    print(report["status"].value_counts())
    print("Report saved.")

if __name__ == "__main__":
    main()
