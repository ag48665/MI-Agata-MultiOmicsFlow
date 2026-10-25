
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


BASE = Path("month3_transcriptomics/single_cell")
INPUT = BASE / "outputs/gse176078_annotated_9566.h5ad"
TABLE_DIR = BASE / "outputs/tables"
FIGURE_DIR = BASE / "outputs/figures"
REPORT_DIR = BASE / "reports"

for folder in (TABLE_DIR, FIGURE_DIR, REPORT_DIR):
    folder.mkdir(parents=True, exist_ok=True)


def main():
    print("=== SAMPLE-LEVEL CELL COMPOSITION ===")

    adata = sc.read_h5ad(INPUT)

    required = ["orig.ident", "subtype", "celltype_annotated"]
    missing = [col for col in required if col not in adata.obs.columns]

    if missing:
        raise ValueError(f"Missing metadata columns: {missing}")

    df = adata.obs[required].copy()
    df = df.dropna(subset=required)

    sample_subtypes = df.groupby(
        "orig.ident", observed=True
    )["subtype"].nunique()

    if (sample_subtypes > 1).any():
        raise ValueError(
            "Some samples have more than one subtype."
        )

    counts = (
        df.groupby(
            ["orig.ident", "subtype", "celltype_annotated"],
            observed=True
        )
        .size()
        .rename("cell_count")
        .reset_index()
    )

    sample_totals = (
        df.groupby("orig.ident", observed=True)
        .size()
        .rename("sample_total")
        .reset_index()
    )

    counts = counts.merge(
        sample_totals,
        on="orig.ident",
        validate="many_to_one"
    )

    counts["percentage"] = (
        counts["cell_count"] / counts["sample_total"] * 100
    )

    # Include zero percentages for cell types absent from a sample.
    samples = df[["orig.ident", "subtype"]].drop_duplicates()
    cell_types = sorted(df["celltype_annotated"].unique())

    complete = samples.merge(
        pd.DataFrame({"celltype_annotated": cell_types}),
        how="cross"
    )

    counts = complete.merge(
        counts,
        on=["orig.ident", "subtype", "celltype_annotated"],
        how="left"
    )

    counts["cell_count"] = counts["cell_count"].fillna(0).astype(int)
    counts["sample_total"] = counts["sample_total"].fillna(
        counts["orig.ident"].map(
            sample_totals.set_index("orig.ident")["sample_total"]
        )
    ).astype(int)
    counts["percentage"] = counts["percentage"].fillna(0.0)

    counts = counts.sort_values(
        ["subtype", "orig.ident", "celltype_annotated"]
    )

    table_path = TABLE_DIR / "celltype_composition_per_sample.csv"
    counts.to_csv(table_path, index=False)

    # One observation per sample for each cell type.
    plot_data = counts.copy()
    subtypes = ["ER+", "HER2+", "TNBC"]
    available = [s for s in subtypes if s in plot_data["subtype"].unique()]

    fig, axes = plt.subplots(
        len(cell_types),
        1,
        figsize=(10, max(3 * len(cell_types), 6)),
        sharex=True
    )

    if len(cell_types) == 1:
        axes = [axes]

    for ax, cell_type in zip(axes, cell_types):
        subset = plot_data[
            plot_data["celltype_annotated"] == cell_type
        ]

        values = [
            subset.loc[
                subset["subtype"] == subtype, "percentage"
            ].to_numpy()
            for subtype in available
        ]

        ax.boxplot(values, tick_labels=available, showfliers=False)

        for index, group_values in enumerate(values, start=1):
            if len(group_values):
                # Deterministic horizontal spacing for visible points.
                positions = [
                    index + (j % 5 - 2) * 0.035
                    for j in range(len(group_values))
                ]
                ax.scatter(
                    positions,
                    group_values,
                    s=22,
                    alpha=0.7,
                    zorder=3
                )

        ax.set_title(str(cell_type))
        ax.set_ylabel("Cells (%)")
        ax.grid(axis="y", alpha=0.25)

    axes[-1].set_xlabel("Breast cancer subtype")
    fig.suptitle(
        "Cell-type composition across individual samples",
        fontsize=14
    )
    fig.tight_layout()

    figure_path = FIGURE_DIR / "celltype_composition_per_sample.png"
    fig.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    report = {
        "dataset": "GSE176078",
        "cells": int(len(df)),
        "samples": int(df["orig.ident"].nunique()),
        "samples_per_subtype": {
            str(k): int(v)
            for k, v in df.groupby(
                "subtype", observed=True
            )["orig.ident"].nunique().items()
        },
        "cell_types": [str(x) for x in cell_types],
        "analysis_level": "sample",
        "statistical_tests_performed": False,
        "note": (
            "Sample-level descriptive analysis. "
            "Biological independence requires source metadata verification."
        )
    }

    report_path = REPORT_DIR / "sample_composition_report.json"
    report_path.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8"
    )

    print(f"Cells: {len(df)}")
    print(f"Samples: {df['orig.ident'].nunique()}")
    print(f"Cell types: {len(cell_types)}")
    print(f"Table: {table_path}")
    print(f"Figure: {figure_path}")
    print(f"Report: {report_path}")
    print("=== SAMPLE COMPOSITION FINISHED ===")


if __name__ == "__main__":
    main()


