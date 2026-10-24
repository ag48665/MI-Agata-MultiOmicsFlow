
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 1. PROJECT PATHS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "figures"

FIGURES.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. ENRICHMENT LIBRARIES
# ============================================================

LIBRARIES = {
    "GO_Biological_Process_2023": {
        "title": "TCGA-BRCA | GO Biological Process",
        "filename": "tcga_brca_go_enrichment.png",
        "color": "#2878B5",
    },
    "KEGG_2021_Human": {
        "title": "TCGA-BRCA | KEGG Pathway Enrichment",
        "filename": "tcga_brca_kegg_enrichment.png",
        "color": "#C44E52",
    },
}


# ============================================================
# 3. VISUALIZATION FUNCTION
# ============================================================

def visualize(library, config):

    input_path = TABLES / f"{library}_enrichment.csv"

    print(f"\nAnalyzing: {library}")

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_path}"
        )

    df = pd.read_csv(input_path)

    required = {
        "Term",
        "Overlap",
        "Adjusted P-value"
    }

    if not required.issubset(df.columns):
        raise ValueError(
            f"Missing required columns in {input_path}"
        )

    # Convert values to numeric
    df["Adjusted P-value"] = pd.to_numeric(
        df["Adjusted P-value"],
        errors="coerce"
    )

    df["Overlap"] = pd.to_numeric(
        df["Overlap"],
        errors="coerce"
    )

    # Remove missing values
    df = df.dropna(
        subset=[
            "Term",
            "Overlap",
            "Adjusted P-value"
        ]
    ).copy()

    # Keep significant pathways
    df = df[
        (df["Adjusted P-value"] < 0.05)
        & (df["Adjusted P-value"] >= 0)
        & (df["Overlap"] > 0)
    ].copy()

    print("Significant pathways:", len(df))

    if df.empty:
        print("No significant pathways found.")
        return 0

    # Select top 20 pathways
    df = df.sort_values(
        "Adjusted P-value"
    ).head(20).copy()

    # Calculate -log10(FDR)
    df["minus_log10_fdr"] = -np.log10(
        df["Adjusted P-value"].clip(
            lower=np.nextafter(0, 1)
        )
    )

    # Remove GO identifiers from labels
    df["Label"] = df["Term"].str.replace(
        r"\s*\(GO:\d+\)$",
        "",
        regex=True
    )

    # Reverse order for horizontal plot
    df = df.iloc[::-1]

    # ========================================================
    # 4. CREATE FIGURE
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(13, 9)
    )

    bars = ax.barh(
        df["Label"],
        df["minus_log10_fdr"],
        color=config["color"],
        alpha=0.85,
        edgecolor="white"
    )

    ax.set_xlabel(
        "-log10(adjusted p-value)",
        fontsize=12
    )

    ax.set_ylabel("")

    ax.set_title(
        config["title"],
        fontsize=16,
        fontweight="bold",
        pad=20
    )

    ax.grid(
        axis="x",
        linestyle="--",
        alpha=0.25
    )

    ax.set_axisbelow(True)

    # ========================================================
    # 5. ADD GENE COUNTS
    # ========================================================

    max_value = df["minus_log10_fdr"].max()

    ax.set_xlim(
        0,
        max_value * 1.18
    )

    for bar, overlap in zip(
        bars,
        df["Overlap"]
    ):

        ax.text(
            bar.get_width() + max_value * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"n={int(overlap)}",
            va="center",
            fontsize=9
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ========================================================
    # 6. SAVE FIGURE
    # ========================================================

    fig.tight_layout()

    output_path = FIGURES / config["filename"]

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("Pathways displayed:", len(df))
    print("Saved:", output_path)

    return len(df)


# ============================================================
# 7. MAIN FUNCTION
# ============================================================

def main():

    print("=" * 50)
    print("TCGA-BRCA PATHWAY VISUALIZATION")
    print("=" * 50)

    total_plots = 0

    for library, config in LIBRARIES.items():

        count = visualize(
            library,
            config
        )

        if count > 0:
            total_plots += 1

    print("\n" + "=" * 50)
    print("PATHWAY VISUALIZATION SUMMARY")
    print("=" * 50)

    print("Libraries analyzed:", len(LIBRARIES))
    print("Plots generated:", total_plots)

    if total_plots == len(LIBRARIES):
        print("STATUS: PASS")
    else:
        print("STATUS: INCOMPLETE")

    print("PATHWAY VISUALIZATION COMPLETED")


# ============================================================
# 8. RUN SCRIPT
# ============================================================

if __name__ == "__main__":
    main()
