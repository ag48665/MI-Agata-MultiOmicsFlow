from pathlib import Path
import json

import scanpy as sc
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Project directories
BASE = Path(__file__).resolve().parents[1]
INPUT = BASE / "outputs" / "gse176078_subset_10000.h5ad"
OUTPUTS = BASE / "outputs"
FIGURES = OUTPUTS / "figures"
REPORTS = BASE / "reports"

FIGURES.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

print("=== GSE176078 SINGLE-CELL QC ===")

# Load data
adata = sc.read_h5ad(INPUT)

print("Cells before QC:", adata.n_obs)
print("Genes before QC:", adata.n_vars)

# Convert existing metadata to numeric values
for column in ["nCount_RNA", "nFeature_RNA", "percent.mito"]:
    if column not in adata.obs.columns:
        raise ValueError(f"Missing metadata column: {column}")

    adata.obs[column] = pd.to_numeric(
        adata.obs[column],
        errors="coerce"
    )

# Save descriptive QC statistics
qc_columns = [
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mito"
]

summary = adata.obs[qc_columns].describe()
summary.to_csv(REPORTS / "qc_summary_before_filtering.csv")

print("\nQC summary:")
print(summary)

# Visualize QC metrics
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, column in zip(axes, qc_columns):
    values = adata.obs[column].dropna()
    ax.hist(values, bins=50)
    ax.set_title(column)
    ax.set_xlabel(column)
    ax.set_ylabel("Number of cells")

plt.tight_layout()
plt.savefig(
    FIGURES / "qc_metrics_before_filtering.png",
    dpi=200
)
plt.close()

# Initial exploratory QC thresholds
min_genes = 200
max_genes = 7500
max_mito = 20.0

# percent.mito in this dataset is expected to be a percentage (0-100).
# Verify its distribution before accepting the thresholds.
mito_max_observed = adata.obs["percent.mito"].max()
print("\nMaximum observed percent.mito:", mito_max_observed)

if pd.isna(mito_max_observed):
    raise ValueError("Mitochondrial percentages are missing.")

if mito_max_observed <= 1:
    raise ValueError(
        "percent.mito may be stored as a fraction (0-1), "
        "not a percentage (0-100). Verify units before filtering."
    )

mask = (
    (adata.obs["nFeature_RNA"] >= min_genes)
    & (adata.obs["nFeature_RNA"] <= max_genes)
    & (adata.obs["percent.mito"] <= max_mito)
    & adata.obs[qc_columns].notna().all(axis=1)
)

cells_passing = int(mask.sum())
cells_removed = int(adata.n_obs - cells_passing)

print("\nCells passing exploratory QC:", cells_passing)
print("Cells excluded:", cells_removed)

report = {
    "dataset": "GSE176078",
    "input_cells": int(adata.n_obs),
    "input_genes": int(adata.n_vars),
    "cells_passing_qc": cells_passing,
    "cells_excluded": cells_removed,
    "thresholds": {
        "min_genes": min_genes,
        "max_genes": max_genes,
        "max_mito_percent": max_mito
    },
    "note": (
        "Exploratory QC only. Thresholds require review "
        "before producing a filtered dataset."
    )
}

(REPORTS / "qc_exploratory_report.json").write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print("\nQC report saved.")
print("QC figure saved.")
print("=== QC FINISHED ===")