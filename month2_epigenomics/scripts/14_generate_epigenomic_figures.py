from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
)

FIGURES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "figures"
)

REPORTS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
)

FIGURES_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

REPORTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

INPUT_FILE = (
    PROCESSED_DIR
    / "clinical_genomic_epigenomic_merged.csv"
)

REPORT_FILE = (
    REPORTS_DIR
    / "epigenomic_figures_report.md"
)


# =====================================
# FUNCTIONS
# =====================================

def load_dataset() -> pd.DataFrame:
    """Load the integrated clinical-genomic-epigenomic dataset."""

    print("Loading integrated dataset...")

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{INPUT_FILE}"
        )

    return pd.read_csv(
        INPUT_FILE,
        low_memory=False,
    )


def identify_methylation_columns(
    df: pd.DataFrame,
) -> list[str]:
    """Identify methylation feature columns."""

    return [
        column
        for column in df.columns
        if str(column).startswith("meth_")
    ]


def plot_beta_distribution(
    df: pd.DataFrame,
    methylation_columns: list[str],
) -> None:
    """Create a histogram of methylation beta values."""

    print("Creating beta-value histogram...")

    values = (
        df[methylation_columns]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .to_numpy()
        .ravel()
    )

    values = values[
        ~np.isnan(values)
    ]

    plt.figure(
        figsize=(8, 5)
    )

    plt.hist(
        values,
        bins=50,
    )

    plt.xlabel("Beta value")
    plt.ylabel("Count")

    plt.title(
        "DNA methylation beta-value distribution"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES_DIR
        / "beta_distribution.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def prepare_methylation_matrix(
    df: pd.DataFrame,
    methylation_columns: list[str],
) -> pd.DataFrame:
    """Prepare numeric methylation data for PCA."""

    matrix = df[
        methylation_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    imputer = SimpleImputer(
        strategy="median"
    )

    imputed = imputer.fit_transform(
        matrix
    )

    return pd.DataFrame(
        imputed,
        columns=methylation_columns,
        index=df.index,
    )


def plot_pca(
    matrix: pd.DataFrame,
) -> tuple[float, float]:
    """Create a PCA plot of methylation profiles."""

    print("Creating PCA plot...")

    pca = PCA(
        n_components=2,
        random_state=42,
    )

    components = pca.fit_transform(
        matrix
    )

    explained = (
        pca.explained_variance_ratio_
        * 100
    )

    plt.figure(
        figsize=(8, 6)
    )

    plt.scatter(
        components[:, 0],
        components[:, 1],
        alpha=0.7,
    )

    plt.xlabel(
        f"PC1 ({explained[0]:.2f}% variance)"
    )

    plt.ylabel(
        f"PC2 ({explained[1]:.2f}% variance)"
    )

    plt.title(
        "PCA of DNA methylation profiles"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES_DIR
        / "methylation_pca.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    return (
        float(explained[0]),
        float(explained[1]),
    )


def create_report(
    dataset: pd.DataFrame,
    methylation_columns: list[str],
    pc1_variance: float,
    pc2_variance: float,
) -> str:
    """Create a Markdown report describing generated figures."""

    report_lines = [
        "# Epigenomic Figures Report",
        "",
        "## Dataset",
        "",
        f"- Integrated patients: {dataset.shape[0]:,}",
        (
            "- Methylation features: "
            f"{len(methylation_columns):,}"
        ),
        "",
        "## Generated figures",
        "",
        "- `beta_distribution.png`",
        "- `methylation_pca.png`",
        "",
        "## PCA variance explained",
        "",
        f"- PC1: {pc1_variance:.2f}%",
        f"- PC2: {pc2_variance:.2f}%",
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Generate epigenomic figures."""

    dataset = load_dataset()

    methylation_columns = (
        identify_methylation_columns(
            dataset
        )
    )

    if not methylation_columns:
        raise ValueError(
            "No methylation columns starting "
            "with 'meth_' were found."
        )

    print(
        f"Integrated rows: "
        f"{dataset.shape[0]:,}"
    )

    print(
        f"Methylation features: "
        f"{len(methylation_columns):,}"
    )

    plot_beta_distribution(
        dataset,
        methylation_columns,
    )

    methylation_matrix = (
        prepare_methylation_matrix(
            dataset,
            methylation_columns,
        )
    )

    pc1_variance, pc2_variance = plot_pca(
        methylation_matrix
    )

    report = create_report(
        dataset=dataset,
        methylation_columns=methylation_columns,
        pc1_variance=pc1_variance,
        pc2_variance=pc2_variance,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")

    print(
        "Beta distribution saved to:\n"
        f"{FIGURES_DIR / 'beta_distribution.png'}"
    )

    print(
        "\nPCA figure saved to:\n"
        f"{FIGURES_DIR / 'methylation_pca.png'}"
    )

    print(
        "\nReport saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()