from pathlib import Path
import json

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("month3_transcriptomics/single_cell")

INPUT = BASE / "outputs/gse176078_annotated_9566.h5ad"
TABLES = BASE / "outputs/tables"
FIGURES = BASE / "outputs/figures"
REPORTS = BASE / "reports"

for folder in [TABLES, FIGURES, REPORTS]:
    folder.mkdir(parents=True, exist_ok=True)

print("=== SINGLE-CELL COMPOSITION ANALYSIS ===")

adata = sc.read_h5ad(INPUT)

required = ["subtype", "celltype_annotated"]

for column in required:
    if column not in adata.obs.columns:
        raise ValueError(f"Missing column: {column}")

df = adata.obs[required].copy()

# Liczba komórek według podtypu i typu komórki.
counts = pd.crosstab(
    df["subtype"],
    df["celltype_annotated"]
)

# Udział procentowy w obrębie każdego podtypu.
percentages = counts.div(
    counts.sum(axis=1),
    axis=0
) * 100

counts.to_csv(
    TABLES / "celltype_counts_by_subtype.csv"
)

percentages.round(2).to_csv(
    TABLES / "celltype_percentages_by_subtype.csv"
)

print("\n=== CELL COUNTS ===")
print(counts.to_string())

print("\n=== CELL PERCENTAGES ===")
print(percentages.round(2).to_string())

# Wykres składu komórkowego.
order = ["ER+", "HER2+", "TNBC"]
available = [x for x in order if x in percentages.index]

plot_data = percentages.loc[available]

ax = plot_data.plot(
    kind="bar",
    stacked=True,
    figsize=(12, 7),
    width=0.7,
    colormap="tab20"
)

ax.set_title(
    "GSE176078: Cell-type composition by breast cancer subtype"
)
ax.set_xlabel("Breast cancer subtype")
ax.set_ylabel("Cell proportion (%)")
ax.set_ylim(0, 100)

ax.legend(
    title="Cell type",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.xticks(rotation=0)
plt.tight_layout()

figure_path = FIGURES / "celltype_composition_by_subtype.png"

plt.savefig(
    figure_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

report = {
    "dataset": "GSE176078",
    "total_cells": int(len(df)),
    "subtype_cell_counts": {
        str(k): int(v)
        for k, v in df["subtype"].value_counts().items()
    },
    "cell_types": int(df["celltype_annotated"].nunique()),
    "analysis": "descriptive cell-level composition",
    "limitation": (
        "Cell-level proportions are descriptive; "
        "patient-level replication is required "
        "for inferential subtype comparisons."
    )
}

report_path = REPORTS / "celltype_composition_report.json"

report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print(f"\nFigure: {figure_path}")
print(f"Report: {report_path}")
print("=== COMPOSITION ANALYSIS FINISHED ===")

