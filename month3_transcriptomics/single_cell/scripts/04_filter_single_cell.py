from pathlib import Path
import json

import scanpy as sc
import pandas as pd

# Project directories
BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = BASE_DIR / "outputs"
REPORT_DIR = BASE_DIR / "reports"

INPUT_FILE = OUTPUT_DIR / "gse176078_subset_10000.h5ad"
OUTPUT_FILE = OUTPUT_DIR / "gse176078_filtered_9566.h5ad"

REPORT_DIR.mkdir(parents=True, exist_ok=True)

print("=== SINGLE-CELL FILTERING ===")

adata = sc.read_h5ad(INPUT_FILE)

print(f"Cells before filtering: {adata.n_obs}")
print(f"Genes before filtering: {adata.n_vars}")

genes = pd.to_numeric(
    adata.obs["nFeature_RNA"], errors="coerce"
)
mito = pd.to_numeric(
    adata.obs["percent.mito"], errors="coerce"
)

mask = (
    (genes >= 200)
    & (genes <= 7500)
    & (mito <= 15)
)

adata_filtered = adata[mask.to_numpy(), :].copy()

print(f"Cells after filtering: {adata_filtered.n_obs}")
print(f"Cells removed: {adata.n_obs - adata_filtered.n_obs}")

# Preserve raw counts for downstream analysis
adata_filtered.layers["counts"] = adata_filtered.X.copy()

adata_filtered.write_h5ad(
    OUTPUT_FILE,
    compression="gzip"
)

report = {
    "dataset": "GSE176078",
    "input_cells": int(adata.n_obs),
    "filtered_cells": int(adata_filtered.n_obs),
    "removed_cells": int(adata.n_obs - adata_filtered.n_obs),
    "genes": int(adata_filtered.n_vars),
    "qc_thresholds": {
        "min_genes": 200,
        "max_genes": 7500,
        "max_percent_mito": 15
    },
    "output_file": OUTPUT_FILE.name
}

report_path = REPORT_DIR / "filtering_report.json"

report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print(f"Saved: {OUTPUT_FILE}")
print(f"Report: {report_path}")
print("=== FILTERING FINISHED ===")