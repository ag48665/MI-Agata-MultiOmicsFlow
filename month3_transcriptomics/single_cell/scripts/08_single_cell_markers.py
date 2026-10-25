from pathlib import Path
import json

import scanpy as sc

BASE = Path("month3_transcriptomics/single_cell")
INPUT = BASE / "outputs/gse176078_umap_leiden_9566.h5ad"
TABLES = BASE / "outputs/tables"
REPORTS = BASE / "reports"

TABLES.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

print("=== SINGLE-CELL MARKER ANALYSIS ===")

adata = sc.read_h5ad(INPUT)

print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# Znajdź kolumnę z klastrami Leiden.
candidates = [
    column for column in adata.obs.columns
    if "leiden" in column.lower()
]

if not candidates:
    raise ValueError("No Leiden cluster column found.")

cluster_key = candidates[0]
print(f"Cluster column: {cluster_key}")

# Analiza markerów na danych log-normalizowanych.
# Jeżeli istnieje warstwa raw, sprawdzamy jej zawartość
# przed automatycznym wykorzystaniem.
sc.tl.rank_genes_groups(
    adata,
    groupby=cluster_key,
    method="wilcoxon",
    use_raw=False,
    n_genes=adata.n_vars,
    pts=True,
)

markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

output_csv = TABLES / "leiden_marker_genes.csv"
markers.to_csv(output_csv, index=False)

top_markers = (
    markers.sort_values(
        ["group", "scores"],
        ascending=[True, False],
    )
    .groupby("group", sort=False)
    .head(10)
)

top_csv = TABLES / "leiden_top10_markers.csv"
top_markers.to_csv(top_csv, index=False)

report = {
    "dataset": "GSE176078",
    "cells": int(adata.n_obs),
    "genes": int(adata.n_vars),
    "cluster_column": cluster_key,
    "clusters": int(adata.obs[cluster_key].nunique()),
    "method": "wilcoxon",
    "marker_rows": int(len(markers)),
}

report_path = REPORTS / "marker_analysis_report.json"
report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8",
)

print(f"Clusters: {report['clusters']}")
print(f"Marker rows: {len(markers)}")
print(f"Saved: {output_csv}")
print(f"Top markers: {top_csv}")

print(f"Report: {report_path}")
print("=== MARKER ANALYSIS FINISHED ===")