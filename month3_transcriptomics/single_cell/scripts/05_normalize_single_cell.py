from pathlib import Path
import json

import scanpy as sc
import numpy as np

BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = BASE_DIR / "outputs"
REPORT_DIR = BASE_DIR / "reports"

INPUT_FILE = OUTPUT_DIR / "gse176078_filtered_9566.h5ad"
OUTPUT_FILE = OUTPUT_DIR / "gse176078_normalized_9566.h5ad"

print("=== SINGLE-CELL NORMALIZATION ===")

adata = sc.read_h5ad(INPUT_FILE)

print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# Preserve original counts in the counts layer
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform normalized expression
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat",
    subset=False
)

hvg_count = int(adata.var["highly_variable"].sum())

print(f"Highly variable genes: {hvg_count}")

# Save normalized data
adata.write_h5ad(
    OUTPUT_FILE,
    compression="gzip"
)

report = {
    "dataset": "GSE176078",
    "cells": int(adata.n_obs),
    "genes": int(adata.n_vars),
    "normalization": "total-count normalization",
    "target_sum": 10000,
    "transformation": "log1p",
    "hvg_method": "seurat",
    "highly_variable_genes": hvg_count,
    "output_file": OUTPUT_FILE.name
}

REPORT_DIR.mkdir(parents=True, exist_ok=True)

report_path = REPORT_DIR / "normalization_report.json"
report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print(f"Saved: {OUTPUT_FILE}")
print(f"Report: {report_path}")
print("=== NORMALIZATION FINISHED ===")
