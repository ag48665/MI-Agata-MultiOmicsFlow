from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("month3_transcriptomics/integration")
INPUT = BASE / "outputs/tables/bulk_single_cell_gene_overlap.csv"
OUTPUT = BASE / "outputs/figures/bulk_single_cell_overlap_by_cluster.png"

df = pd.read_csv(INPUT)

df["bulk_direction"] = df["bulk_log2foldchange"].apply(
    lambda x: "UP in Bulk" if x > 0 else (
        "DOWN in Bulk" if x < 0 else "ZERO"
    )
)

counts = pd.crosstab(
    df["leiden_cluster"].astype(str),
    df["bulk_direction"]
)

# Show all 25 Leiden clusters, including those with zero matches.
counts = counts.reindex(
    [str(i) for i in range(25)],
    fill_value=0
)

columns = ["DOWN in Bulk", "UP in Bulk"]
counts = counts.reindex(columns=columns, fill_value=0)

fig, ax = plt.subplots(figsize=(13, 6))

counts.plot(
    kind="bar",
    stacked=True,
    ax=ax,
    color=["#377eb8", "#e41a1c"],
    width=0.8
)

ax.set_title(
    "Bulk RNA-seq and Single-cell RNA-seq gene overlap by Leiden cluster"
)
ax.set_xlabel("Leiden cluster")
ax.set_ylabel("Number of matched cluster-gene pairs")
ax.legend(title="Bulk RNA-seq direction")
ax.grid(axis="y", alpha=0.25)

plt.xticks(rotation=0)
plt.tight_layout()

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close()

print("=== INTEGRATION FIGURE ===")
print("Figure:", OUTPUT)
print("Total overlap rows:", int(counts.to_numpy().sum()))
print("Clusters with matches:", int((counts.sum(axis=1) > 0).sum()))
print("Clusters without matches:", int((counts.sum(axis=1) == 0).sum()))
print("=== FIGURE FINISHED ===")
