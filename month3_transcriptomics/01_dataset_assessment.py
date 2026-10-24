
from pathlib import Path
import pandas as pd

# Project directories
ROOT = Path(__file__).resolve().parents[2]
MONTH3 = ROOT / "month3_transcriptomics"
REPORTS = MONTH3 / "reports"

REPORTS.mkdir(parents=True, exist_ok=True)

# Identify existing data files
extensions = {
    ".csv", ".tsv", ".txt", ".parquet",
    ".h5ad", ".h5", ".rds", ".gz"
}

records = []

for file in ROOT.rglob("*"):
    if not file.is_file():
        continue

    if ".venv" in file.parts or ".git" in file.parts:
        continue

    if file.suffix.lower() in extensions:
        records.append({
            "file": str(file.relative_to(ROOT)),
            "size_mb": round(
                file.stat().st_size / (1024 ** 2), 3
            ),
            "extension": file.suffix.lower()
        })

df = pd.DataFrame(
    records,
    columns=["file", "size_mb", "extension"]
)

output = REPORTS / "local_dataset_inventory.csv"
df.to_csv(output, index=False)

print("TCGA-BRCA Month 3 - Dataset Assessment")
print("-------------------------------------")
print(f"Files found: {len(df)}")
print(f"Report saved: {output}")

if not df.empty:
    print(df.to_string(index=False))
