
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    TABLES / "tcga_brca_differential_expression.csv"
)

required = {"log2FoldChange", "padj"}
if not required.issubset(df.columns):
    raise ValueError("Missing required DEG columns.")

df = df.copy()

# Avoid log10(0) for extremely small adjusted p-values
positive = df.loc[df["padj"] > 0, "padj"]

if positive.empty:
    raise ValueError("No positive adjusted p-values.")

minimum = max(
    positive.min() / 10,
    np.finfo(float).tiny
)

df["neg_log10_padj"] = -np.log10(
    df["padj"].clip(lower=minimum)
)

up = (
    (df["padj"] < 0.05) &
    (df["log2FoldChange"] > 1)
)

down = (
    (df["padj"] < 0.05) &
    (df["log2FoldChange"] < -1)
)

other = ~(up | down)

fig, ax = plt.subplots(figsize=(11, 7))

ax.scatter(
    df.loc[other, "log2FoldChange"],
    df.loc[other, "neg_log10_padj"],
    color="lightgray", s=5, alpha=0.5,
    label=f"Other ({other.sum():,})"
)

ax.scatter(
    df.loc[down, "log2FoldChange"],
    df.loc[down, "neg_log10_padj"],
    color="royalblue", s=7, alpha=0.65,
    label=f"Downregulated ({down.sum():,})"
)

ax.scatter(
    df.loc[up, "log2FoldChange"],
    df.loc[up, "neg_log10_padj"],
    color="firebrick", s=7, alpha=0.65,
    label=f"Upregulated ({up.sum():,})"
)

ax.axvline(-1, linestyle="--", color="black", linewidth=0.8)
ax.axvline(1, linestyle="--", color="black", linewidth=0.8)
ax.axhline(
    -np.log10(0.05),
    linestyle="--", color="black", linewidth=0.8
)

ax.set_xlabel("log2 Fold Change (Tumor / Normal)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("TCGA-BRCA | Paired Differential Expression")
ax.legend()
ax.grid(alpha=0.15)

fig.tight_layout()

output = FIGURES / "tcga_brca_volcano_plot.png"
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("VOLCANO PLOT")
print("Genes:", len(df))
print("Upregulated:", up.sum())
print("Downregulated:", down.sum())
print("Saved:", output)
print("STATUS: PASS")
