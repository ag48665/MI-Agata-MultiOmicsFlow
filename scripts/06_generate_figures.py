from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# =====================================
# ŚCIEŻKI
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"
FIGURES_DIR = PROJECT_ROOT / "outputs" / "figures"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

MERGED_FILE = PROCESSED_DIR / "clinical_genomics_merged.csv"
TOP_GENES_FILE = TABLES_DIR / "top_mutated_genes.csv"
SUBTYPE_FILE = TABLES_DIR / "cancer_subtype_summary.csv"
VARIANT_FILE = TABLES_DIR / "variant_classification_summary.csv"


# =====================================
# FUNKCJE POMOCNICZE
# =====================================

def load_csv(file_path: Path) -> pd.DataFrame:
    """Wczytuje plik CSV i sprawdza, czy istnieje."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{file_path}"
        )

    return pd.read_csv(file_path, low_memory=False)


def save_figure(file_name: str) -> None:
    """Zapisuje aktualny wykres i zamyka figurę."""

    output_path = FIGURES_DIR / file_name

    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print(f"Saved: {output_path}")


def plot_age_distribution(merged: pd.DataFrame) -> None:
    """Tworzy histogram wieku pacjentów."""

    if "age" not in merged.columns:
        print("Skipped age plot: column 'age' not found.")
        return

    age = pd.to_numeric(
        merged["age"],
        errors="coerce"
    ).dropna()

    if age.empty:
        print("Skipped age plot: no valid age values.")
        return

    plt.figure(figsize=(9, 6))
    plt.hist(age, bins=20, edgecolor="black")
    plt.xlabel("Age")
    plt.ylabel("Number of patients")
    plt.title("Age distribution in TCGA-BRCA")

    save_figure("age_distribution.png")


def plot_mutations_per_patient(
    merged: pd.DataFrame
) -> None:
    """Tworzy histogram liczby mutacji na pacjenta."""

    if "mutation_count" not in merged.columns:
        print(
            "Skipped mutation count plot: "
            "column 'mutation_count' not found."
        )
        return

    mutation_count = pd.to_numeric(
        merged["mutation_count"],
        errors="coerce"
    ).dropna()

    if mutation_count.empty:
        print(
            "Skipped mutation count plot: "
            "no valid values."
        )
        return

    upper_limit = mutation_count.quantile(0.99)
    displayed_values = mutation_count[
        mutation_count <= upper_limit
    ]

    plt.figure(figsize=(9, 6))
    plt.hist(
        displayed_values,
        bins=30,
        edgecolor="black"
    )
    plt.xlabel("Mutation count per patient")
    plt.ylabel("Number of patients")
    plt.title(
        "Mutations per patient in TCGA-BRCA\n"
        "(values up to the 99th percentile)"
    )

    save_figure("mutations_per_patient.png")


def plot_top_mutated_genes(
    top_genes: pd.DataFrame
) -> None:
    """Tworzy wykres najczęściej mutowanych genów."""

    required_columns = {
        "gene",
        "patient_count"
    }

    missing = required_columns.difference(
        top_genes.columns
    )

    if missing:
        print(
            "Skipped top genes plot. "
            f"Missing columns: {sorted(missing)}"
        )
        return

    plot_data = (
        top_genes
        .sort_values(
            "patient_count",
            ascending=True
        )
        .tail(20)
    )

    plt.figure(figsize=(10, 8))
    plt.barh(
        plot_data["gene"],
        plot_data["patient_count"]
    )
    plt.xlabel("Number of patients")
    plt.ylabel("Gene")
    plt.title("Top mutated genes in TCGA-BRCA")

    save_figure("top_mutated_genes_barplot.png")


def plot_subtype_distribution(
    subtype_summary: pd.DataFrame
) -> None:
    """Tworzy wykres rozkładu podtypów nowotworu."""

    subtype_column = None

    for candidate in [
        "subtype",
        "cancer_type_detailed",
        "cancer_type_detailed_sample"
    ]:
        if candidate in subtype_summary.columns:
            subtype_column = candidate
            break

    if subtype_column is None:
        print(
            "Skipped subtype plot: "
            "no subtype column found."
        )
        return

    if "count" not in subtype_summary.columns:
        print(
            "Skipped subtype plot: "
            "column 'count' not found."
        )
        return

    plot_data = (
        subtype_summary
        .sort_values("count", ascending=True)
    )

    plt.figure(figsize=(10, 7))
    plt.barh(
        plot_data[subtype_column].astype(str),
        plot_data["count"]
    )
    plt.xlabel("Number of patients")
    plt.ylabel("Subtype")
    plt.title("Breast cancer subtype distribution")

    save_figure("cancer_subtype_distribution.png")


def plot_variant_classification(
    variant_summary: pd.DataFrame
) -> None:
    """Tworzy wykres typów wariantów."""

    required_columns = {
        "variant_classification",
        "count"
    }

    missing = required_columns.difference(
        variant_summary.columns
    )

    if missing:
        print(
            "Skipped variant plot. "
            f"Missing columns: {sorted(missing)}"
        )
        return

    plot_data = (
        variant_summary
        .head(15)
        .sort_values("count", ascending=True)
    )

    plt.figure(figsize=(11, 8))
    plt.barh(
        plot_data["variant_classification"].astype(str),
        plot_data["count"]
    )
    plt.xlabel("Number of mutation records")
    plt.ylabel("Variant classification")
    plt.title(
        "Distribution of mutation classifications"
    )

    save_figure(
        "variant_classification_distribution.png"
    )


def plot_mutation_burden_by_subtype(
    merged: pd.DataFrame
) -> None:
    """Tworzy boxplot mutacji według podtypu raka."""

    if "mutation_count" not in merged.columns:
        print(
            "Skipped burden-by-subtype plot: "
            "mutation_count not found."
        )
        return

    subtype_column = None

    for candidate in [
        "subtype",
        "subtype_patient",
        "cancer_type_detailed",
        "cancer_type_detailed_sample"
    ]:
        if candidate in merged.columns:
            subtype_column = candidate
            break

    if subtype_column is None:
        print(
            "Skipped burden-by-subtype plot: "
            "subtype column not found."
        )
        return

    plot_data = merged[
        [subtype_column, "mutation_count"]
    ].copy()

    plot_data["mutation_count"] = pd.to_numeric(
        plot_data["mutation_count"],
        errors="coerce"
    )

    plot_data = plot_data.dropna()

    if plot_data.empty:
        print(
            "Skipped burden-by-subtype plot: "
            "no valid data."
        )
        return

    groups = []
    labels = []

    for subtype, group in plot_data.groupby(
        subtype_column
    ):
        if len(group) >= 5:
            groups.append(
                group["mutation_count"].values
            )
            labels.append(str(subtype))

    if not groups:
        print(
            "Skipped burden-by-subtype plot: "
            "not enough observations."
        )
        return

    plt.figure(figsize=(11, 7))
    plt.boxplot(
    groups,
    tick_labels=labels,
    showfliers=False
  )
    plt.xlabel("Subtype")
    plt.ylabel("Mutation count")
    plt.title(
        "Mutation burden by breast cancer subtype"
    )
    plt.xticks(rotation=45, ha="right")

    save_figure(
        "mutation_burden_by_clinical_group.png"
    )


# =====================================
# MAIN
# =====================================

def main() -> None:
    print("Loading input data...")

    merged = load_csv(MERGED_FILE)
    top_genes = load_csv(TOP_GENES_FILE)
    subtype_summary = load_csv(SUBTYPE_FILE)
    variant_summary = load_csv(VARIANT_FILE)

    print("Generating figures...")

    plot_age_distribution(merged)
    plot_mutations_per_patient(merged)
    plot_top_mutated_genes(top_genes)
    plot_subtype_distribution(subtype_summary)
    plot_variant_classification(variant_summary)
    plot_mutation_burden_by_subtype(merged)

    print("\nFinished successfully.")
    print(f"Figures directory: {FIGURES_DIR}")


if __name__ == "__main__":
    main()