from pathlib import Path
import json

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("month3_transcriptomics/single_cell")

INPUT = BASE / "outputs/gse176078_umap_leiden_9566.h5ad"
OUTPUT = BASE / "outputs/gse176078_annotated_9566.h5ad"

FIGURES = BASE / "outputs/figures"
TABLES = BASE / "outputs/tables"
REPORTS = BASE / "reports"

for folder in [FIGURES, TABLES, REPORTS]:
    folder.mkdir(parents=True, exist_ok=True)

print("=== SINGLE-CELL ANNOTATION ===")

adata = sc.read_h5ad(INPUT)

required = ["leiden", "celltype_major"]

for column in required:
    if column not in adata.obs.columns:
        raise ValueError(f"Missing column: {column}")

if "X_umap" not in adata.obsm:
    raise ValueError("UMAP coordinates are missing.")

# Zachowujemy oryginalne adnotacje na poziomie komórek.
adata.obs["celltype_annotated"] = (
    adata.obs["celltype_major"]
    .astype("category")
)

# Porównanie klastrów Leiden z adnotacjami.
table = pd.crosstab(
    adata.obs["leiden"],
    adata.obs["celltype_annotated"]
)

table.to_csv(TABLES / "leiden_celltype_crosstab.csv")

summary = pd.DataFrame({
    "cells": table.sum(axis=1),
    "dominant_celltype": table.idxmax(axis=1),
    "dominant_cells": table.max(axis=1),
})

summary["dominant_pct"] = (
    100 * summary["dominant_cells"] / summary["cells"]
).round(2)

summary.to_csv(
    TABLES / "leiden_celltype_annotation_summary.csv"
)

print("\nAnnotation summary:")
print(summary.to_string())

# UMAP z typami komórek.
sc.pl.umap(
    adata,
    color="celltype_annotated",
    legend_loc="right margin",
    frameon=False,
    show=False,
)

plt.savefig(
    FIGURES / "umap_celltype_annotated.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close("all")

adata.write_h5ad(OUTPUT)

report = {
    "dataset": "GSE176078",
    "cells": int(adata.n_obs),
    "genes": int(adata.n_vars),
    "clusters": int(adata.obs["leiden"].nunique()),
    "cell_types": int(
        adata.obs["celltype_annotated"].nunique()
    ),
    "annotation_source": "GSE176078 celltype_major",
    "method": "original cell-level annotations",
    "clusters_below_90pct_dominance": (
        summary.index[
            summary["dominant_pct"] < 90
        ].astype(str).tolist()
    ),
}

report_path = REPORTS / "annotation_report.json"

report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8",
)

print(f"\nCells: {adata.n_obs}")
print(f"Clusters: {report['clusters']}")
print(f"Cell types: {report['cell_types']}")
print(f"Saved: {OUTPUT}")
print(f"Figure: {FIGURES / 'umap_celltype_annotated.png'}")
print(f"Report: {report_path}")
print("=== ANNOTATION FINISHED ===")

