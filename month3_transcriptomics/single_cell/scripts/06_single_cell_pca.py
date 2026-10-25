from pathlib import Path
import json

import scanpy as sc
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = BASE_DIR / "outputs"
REPORT_DIR = BASE_DIR / "reports"
FIGURE_DIR = OUTPUT_DIR / "figures"

INPUT_FILE = OUTPUT_DIR / "gse176078_normalized_9566.h5ad"
OUTPUT_FILE = OUTPUT_DIR / "gse176078_pca_9566.h5ad"

REPORT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

print("=== SINGLE-CELL PCA ===")

adata = sc.read_h5ad(INPUT_FILE)

print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

hvg_count = int(adata.var["highly_variable"].sum())
print(f"Highly variable genes: {hvg_count}")

if hvg_count < 50:
    raise ValueError("Too few highly variable genes for PCA.")

# Work only with highly variable genes for PCA.
# Keep the full normalized dataset unchanged.
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# Scale expression and limit extreme values
sc.pp.scale(adata_hvg, max_value=10)

# Principal component analysis
sc.tl.pca(
    adata_hvg,
    n_comps=50,
    svd_solver="arpack",
    random_state=42
)

# Store PCA coordinates in the full AnnData object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

variance_ratio = adata_hvg.uns["pca"]["variance_ratio"]

adata.uns["pca"] = {
    "variance": adata_hvg.uns["pca"]["variance"].copy(),
    "variance_ratio": variance_ratio.copy(),
    "params": dict(adata_hvg.uns["pca"]["params"])
}

print(f"Principal components: {adata.obsm['X_pca'].shape[1]}")
print(
    "Variance explained by first 10 PCs: "
    f"{np.sum(variance_ratio[:10]) * 100:.2f}%"
)

# Save explained-variance plot
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(
    np.arange(1, len(variance_ratio) + 1),
    variance_ratio * 100,
    marker="o",
    markersize=3
)

ax.set_xlabel("Principal component")
ax.set_ylabel("Explained variance (%)")
ax.set_title("GSE176078 - PCA explained variance")
ax.grid(alpha=0.3)

fig.tight_layout()

figure_path = FIGURE_DIR / "pca_explained_variance.png"
fig.savefig(figure_path, dpi=200)
plt.close(fig)

# Save PCA dataset
adata.write_h5ad(OUTPUT_FILE, compression="gzip")

report = {
    "dataset": "GSE176078",
    "cells": int(adata.n_obs),
    "genes": int(adata.n_vars),
    "highly_variable_genes": hvg_count,
    "principal_components": 50,
    "variance_explained_first_10_pct": round(
        float(np.sum(variance_ratio[:10]) * 100), 2
    ),
    "random_state": 42,
    "output_file": OUTPUT_FILE.name
}

report_path = REPORT_DIR / "pca_report.json"

report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print(f"Saved: {OUTPUT_FILE}")
print(f"Figure: {figure_path}")
print(f"Report: {report_path}")
print("=== PCA FINISHED ===")

