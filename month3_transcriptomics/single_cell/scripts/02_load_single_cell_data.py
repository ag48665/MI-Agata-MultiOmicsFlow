
from pathlib import Path
from importlib.metadata import version
import gzip
import json
import random

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# MI-Agata-MultiOmicsFlow
# Month 3 - Single-cell RNA-seq
# Dataset: GEO GSE176078

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data" / "Wu_etal_2021_BRCA_scRNASeq"
OUTPUT = BASE / "outputs"
REPORTS = BASE / "reports"

N_CELLS = 10000
SEED = 42

OUTPUT.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

mtx_path = DATA / "count_matrix_sparse.mtx"
genes_path = DATA / "count_matrix_genes.tsv"
barcodes_path = DATA / "count_matrix_barcodes.tsv"
metadata_path = DATA / "metadata.csv"

print("=== GSE176078 SINGLE-CELL LOADING ===", flush=True)

for path in [mtx_path, genes_path, barcodes_path, metadata_path]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")

genes = pd.read_csv(
    genes_path, sep="\t", header=None, dtype=str
).iloc[:, 0].tolist()

barcodes = pd.read_csv(
    barcodes_path, sep="\t", header=None, dtype=str
).iloc[:, 0].tolist()

metadata = pd.read_csv(
    metadata_path, index_col=0, low_memory=False
)
metadata.index = metadata.index.astype(str)

if len(set(barcodes)) != len(barcodes):
    raise ValueError("Duplicate barcodes detected")

if not metadata.index.is_unique:
    raise ValueError("Duplicate metadata identifiers detected")

missing = set(barcodes) - set(metadata.index)
if missing:
    raise ValueError(
        f"Metadata missing for {len(missing)} cells"
    )

n_genes = len(genes)
n_cells = len(barcodes)

print(f"Genes: {n_genes}", flush=True)
print(f"Cells: {n_cells}", flush=True)

rng = random.Random(SEED)
selected = sorted(
    rng.sample(range(n_cells), min(N_CELLS, n_cells))
)

lookup = {
    original: new
    for new, original in enumerate(selected)
}

selected_barcodes = [barcodes[i] for i in selected]

# Matrix Market: rows = genes, columns = cells.
# Read line by line to avoid loading the full matrix.

rows = []
cols = []
values = []

with mtx_path.open("rt", encoding="utf-8") as f:
    header = f.readline().strip()

    if not header.startswith("%%MatrixMarket matrix coordinate"):
        raise ValueError("Unexpected Matrix Market format")

    if "integer" not in header.lower():
        raise ValueError("Expected integer count matrix")

    line = f.readline()
    while line.startswith("%"):
        line = f.readline()

    mtx_genes, mtx_cells, expected_nnz = map(
        int, line.split()
    )

    if (mtx_genes, mtx_cells) != (n_genes, n_cells):
        raise ValueError("Matrix dimensions do not match")

    processed = 0

    for line in f:
        parts = line.split()
        if not parts:
            continue

        gene_id = int(parts[0]) - 1
        cell_id = int(parts[1]) - 1
        count = int(parts[2])

        processed += 1

        new_cell_id = lookup.get(cell_id)
        if new_cell_id is not None:
            rows.append(new_cell_id)
            cols.append(gene_id)
            values.append(count)

        if processed % 20000000 == 0:
            print(
                f"Processed {processed:,} entries",
                flush=True
            )

    if processed != expected_nnz:
        raise ValueError(
            f"Expected {expected_nnz} entries, got {processed}"
        )

print("Building sparse matrix...", flush=True)

X = sp.coo_matrix(
    (
        np.asarray(values, dtype=np.int32),
        (
            np.asarray(rows, dtype=np.int32),
            np.asarray(cols, dtype=np.int32),
        ),
    ),
    shape=(len(selected), n_genes),
    dtype=np.int32,
).tocsr()

X.sum_duplicates()

obs = metadata.loc[selected_barcodes].copy()
var = pd.DataFrame(index=pd.Index(genes, name="gene"))

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.var_names_make_unique()

adata.uns["dataset_id"] = "GSE176078"
adata.uns["random_seed"] = SEED
adata.uns["sampling_method"] = "random_without_replacement"

output_path = OUTPUT / "gse176078_subset_10000.h5ad"
adata.write_h5ad(output_path, compression="gzip")

report = {
    "dataset": "GSE176078",
    "total_cells": n_cells,
    "selected_cells": adata.n_obs,
    "genes": adata.n_vars,
    "nonzero_entries": int(adata.X.nnz),
    "seed": SEED,
    "anndata_version": version("anndata"),
    "output": str(output_path),
}

report_path = REPORTS / "subset_loading_report.json"
report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print("=== LOADING FINISHED ===", flush=True)
print(f"Shape: {adata.shape}", flush=True)
print(f"Output: {output_path}", flush=True)
print(f"Report: {report_path}", flush=True)

