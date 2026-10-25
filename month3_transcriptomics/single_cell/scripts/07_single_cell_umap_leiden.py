from pathlib import Path
from datetime import datetime, timezone
import json

import scanpy as sc

# ==========================================
# MI-Agata-MultiOmicsFlow
# Month 3 - Single-cell RNA-seq
# Step 07: UMAP and Leiden clustering
# Dataset: GEO GSE176078
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[1]

INPUT_FILE = BASE_DIR / "outputs" / "gse176078_pca_9566.h5ad"
OUTPUT_FILE = BASE_DIR / "outputs" / "gse176078_umap_leiden_9566.h5ad"

FIGURE_DIR = BASE_DIR / "outputs" / "figures"
REPORT_DIR = BASE_DIR / "reports"

N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.5
RANDOM_STATE = 42


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n=== GSE176078 UMAP AND LEIDEN ===")

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Missing PCA file: {INPUT_FILE}")

    adata = sc.read_h5ad(INPUT_FILE)

    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")

    if "X_pca" not in adata.obsm:
        raise ValueError("PCA coordinates are missing: X_pca")

    available_pcs = adata.obsm["X_pca"].shape[1]

    if available_pcs < N_PCS:
        raise ValueError(
            f"Requested {N_PCS} PCs, but only {available_pcs} available"
        )

    print(f"Using PCA components: {N_PCS}")
    print(f"Neighbors: {N_NEIGHBORS}")

    sc.pp.neighbors(
        adata,
        n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS,
        use_rep="X_pca",
        random_state=RANDOM_STATE,
    )

    print("Computing UMAP...")

    sc.tl.umap(
        adata,
        random_state=RANDOM_STATE,
    )

    print("Computing Leiden clusters...")

    sc.tl.leiden(
        adata,
        resolution=LEIDEN_RESOLUTION,
        key_added="leiden",
        flavor="igraph",
        directed=False,
        random_state=RANDOM_STATE,
    )

    cluster_counts = adata.obs["leiden"].value_counts().sort_index()

    print(f"Clusters found: {len(cluster_counts)}")
    print(cluster_counts.to_string())

    adata.write_h5ad(OUTPUT_FILE)

    figure_path = FIGURE_DIR / "umap_leiden.png"

    sc.pl.umap(
        adata,
        color="leiden",
        legend_loc="on data",
        show=False,
        return_fig=True,
    ).savefig(
        figure_path,
        dpi=200,
        bbox_inches="tight",
    )

    report = {
        "dataset": "GSE176078",
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
        "leiden_resolution": LEIDEN_RESOLUTION,
        "random_state": RANDOM_STATE,
        "n_clusters": int(len(cluster_counts)),
        "cluster_sizes": {
            str(k): int(v)
            for k, v in cluster_counts.items()
        },
        "input_file": INPUT_FILE.name,
        "output_file": OUTPUT_FILE.name,
    }

    report_path = REPORT_DIR / "umap_leiden_report.json"

    report_path.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Figure: {figure_path}")
    print(f"Report: {report_path}")
    print("=== UMAP AND LEIDEN FINISHED ===")


if __name__ == "__main__":
    main()