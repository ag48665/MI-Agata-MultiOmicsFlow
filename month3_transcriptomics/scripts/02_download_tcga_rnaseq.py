
from pathlib import Path
import requests
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data" / "raw"
OUT.mkdir(parents=True, exist_ok=True)

url = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content": [
        {
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TCGA-BRCA"]
            }
        },
        {
            "op": "=",
            "content": {
                "field": "data_type",
                "value": ["Gene Expression Quantification"]
            }
        },
        {
            "op": "=",
            "content": {
                "field": "analysis.workflow_type",
                "value": ["STAR - Counts"]
            }
        }
    ]
}

fields = ",".join([
    "file_id",
    "file_name",
    "file_size",
    "access",
    "cases.submitter_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type"
])

response = requests.post(
    url,
    json={
        "filters": filters,
        "fields": fields,
        "format": "JSON",
        "size": 2000
    },
    timeout=120
)

response.raise_for_status()
data = response.json()["data"]

records = []

for item in data["hits"]:
    cases = item.get("cases", [])

    if not cases:
        records.append({
            "file_id": item["file_id"],
            "file_name": item["file_name"],
            "file_size_mb": round(
                item.get("file_size", 0) / 1048576, 2
            ),
            "access": item.get("access"),
            "patient_id": None,
            "sample_id": None,
            "sample_type": None
        })

    for case in cases:
        samples = case.get("samples", []) or [{}]

        for sample in samples:
            records.append({
                "file_id": item["file_id"],
                "file_name": item["file_name"],
                "file_size_mb": round(
                    item.get("file_size", 0) / 1048576, 2
                ),
                "access": item.get("access"),
                "patient_id": case.get("submitter_id"),
                "sample_id": sample.get("submitter_id"),
                "sample_type": sample.get("sample_type")
            })

df = pd.DataFrame(records)

output = OUT / "tcga_brca_rnaseq_metadata.csv"
df.to_csv(output, index=False)

print("TCGA-BRCA RNA-seq metadata")
print("-------------------------")
print("Files reported by GDC:", data["pagination"]["total"])
print("Metadata rows:", len(df))
print("Unique file IDs:", df["file_id"].nunique())
print("Unique patients:", df["patient_id"].nunique())
print()
print("Sample types:")
print(df["sample_type"].value_counts(dropna=False))
print()
print("Total file size (GB):",
      round(
          df.drop_duplicates("file_id")["file_size_mb"].sum()
          / 1024, 2
      ))
print()
print("Saved:", output)
