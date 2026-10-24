
from pathlib import Path
import pandas as pd
import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ==========================================
# PATHS
# ==========================================

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
OUTPUT = ROOT / "outputs" / "tables"
OUTPUT.mkdir(parents=True, exist_ok=True)

COUNTS_FILE = PROCESSED / "tcga_brca_raw_counts.csv"
META_FILE = PROCESSED / "tcga_brca_sample_metadata.csv"

# ==========================================
# LOAD DATA
# ==========================================

print("Loading RNA-seq data...")

counts = pd.read_csv(COUNTS_FILE, index_col=0)
metadata = pd.read_csv(META_FILE)

metadata = metadata.set_index("sample_id")

# Convert genes x samples to samples x genes
counts = counts.T

# Match sample order
if not counts.index.is_unique:
    raise ValueError("Duplicate sample IDs in count matrix.")

if not metadata.index.is_unique:
    raise ValueError("Duplicate sample IDs in metadata.")

if set(counts.index) != set(metadata.index):
    raise ValueError("Count matrix and metadata do not match.")

counts = counts.loc[metadata.index]

print("Samples:", counts.shape[0])
print("Genes:", counts.shape[1])

# ==========================================
# VALIDATE PAIRED DESIGN
# ==========================================

metadata["condition"] = metadata["sample_type"].map({
    "Primary Tumor": "Tumor",
    "Solid Tissue Normal": "Normal"
})

if metadata["condition"].isna().any():
    raise ValueError("Unknown sample type.")

pairs = metadata.groupby("patient_id")["condition"].agg(
    lambda x: set(x) == {"Tumor", "Normal"} and len(x) == 2
)

if not pairs.all():
    raise ValueError("Incomplete or duplicated patient pairs.")

metadata["condition"] = pd.Categorical(
    metadata["condition"],
    categories=["Normal", "Tumor"]
)

metadata["patient_id"] = metadata["patient_id"].astype(str)

print("Validated patient pairs:", len(pairs))

# ==========================================
# FILTER LOW-COUNT GENES
# ==========================================

if counts.isna().any().any():
    raise ValueError("Missing counts detected.")

values = counts.to_numpy()

if not np.isfinite(values).all():
    raise ValueError("Non-finite counts detected.")

if (values < 0).any() or (values != np.floor(values)).any():
    raise ValueError("Counts must be nonnegative integers.")

counts = counts.astype(np.int64)

# Keep genes with >= 10 counts in >= 10 samples
keep = (counts >= 10).sum(axis=0) >= 10
counts = counts.loc[:, keep].copy()

print("Genes after filtering:", counts.shape[1])

if counts.shape[1] == 0:
    raise ValueError("No genes remain after filtering.")

# ==========================================
# DIFFERENTIAL EXPRESSION
# ==========================================

print("Running paired differential expression...")

dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design="~ patient_id + condition",
    refit_cooks=True,
    n_cpus=2
)

dds.deseq2()

stats = DeseqStats(
    dds,
    contrast=["condition", "Tumor", "Normal"],
    n_cpus=2
)

stats.summary()

results = stats.results_df.copy()
results.index.name = "gene_id"

# ==========================================
# SAVE RESULTS
# ==========================================

results.to_csv(
    OUTPUT / "tcga_brca_differential_expression.csv"
)

significant = results.loc[
    (results["padj"] < 0.05)
    & (results["log2FoldChange"].abs() > 1)
].copy()

significant.to_csv(
    OUTPUT / "tcga_brca_significant_genes.csv"
)

up = (significant["log2FoldChange"] > 1).sum()
down = (significant["log2FoldChange"] < -1).sum()

print("\nDIFFERENTIAL EXPRESSION RESULTS")
print("Genes tested:", len(results))
print("Significant genes:", len(significant))
print("Upregulated in tumor:", up)
print("Downregulated in tumor:", down)
print("STATUS: PASS")
