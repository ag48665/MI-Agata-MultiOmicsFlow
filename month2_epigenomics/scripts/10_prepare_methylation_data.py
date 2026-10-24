from pathlib import Path

import pandas as pd


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA"
REPORTS_DIR = PROJECT_ROOT / "outputs" / "reports"

REPORTS_DIR.mkdir(parents=True, exist_ok=True)

METHYLATION_450_FILE = (
    RAW_DATA_DIR / "data_methylation_hm450.txt"
)

METHYLATION_27_450_FILE = (
    RAW_DATA_DIR / "data_methylation_hm27_hm450_merged.txt"
)

REPORT_FILE = (
    REPORTS_DIR / "methylation_data_overview.md"
)


# =====================================
# FUNCTIONS
# =====================================

def check_file(file_path: Path) -> None:
    """Check whether the required file exists."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{file_path}"
        )


def load_methylation_data(
    file_path: Path
) -> pd.DataFrame:
    """Load TCGA methylation data."""

    return pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        low_memory=False,
    )


def identify_sample_columns(
    methylation: pd.DataFrame
) -> list[str]:
    """
    Identify TCGA sample columns.

    Sample columns are expected to start with 'TCGA-'.
    """

    return [
        column
        for column in methylation.columns
        if str(column).upper().startswith("TCGA-")
    ]


def create_methylation_summary(
    methylation: pd.DataFrame,
    file_path: Path
) -> str:
    """Create a Markdown summary of the methylation dataset."""

    sample_columns = identify_sample_columns(
        methylation
    )

    metadata_columns = [
        column
        for column in methylation.columns
        if column not in sample_columns
    ]

    if sample_columns:
        sample_data = methylation[
            sample_columns
        ].apply(
            pd.to_numeric,
            errors="coerce",
        )

        missing_values = int(
            sample_data.isna().sum().sum()
        )

        total_values = int(
            sample_data.shape[0]
            * sample_data.shape[1]
        )

        missing_percent = (
            missing_values / total_values * 100
            if total_values
            else 0
        )

        minimum_value = sample_data.min().min()
        maximum_value = sample_data.max().max()
    else:
        missing_values = 0
        total_values = 0
        missing_percent = 0
        minimum_value = pd.NA
        maximum_value = pd.NA

    summary = [
        "# TCGA-BRCA Methylation Data Overview",
        "",
        "## Input file",
        "",
        f"- File name: `{file_path.name}`",
        f"- File path: `{file_path}`",
        "",
        "## Dataset dimensions",
        "",
        f"- Number of rows/features: {methylation.shape[0]:,}",
        f"- Total number of columns: {methylation.shape[1]:,}",
        f"- Number of sample columns: {len(sample_columns):,}",
        f"- Number of metadata columns: {len(metadata_columns):,}",
        "",
        "## Metadata columns",
        "",
        "```text",
        "\n".join(
            str(column)
            for column in metadata_columns
        ),
        "```",
        "",
        "## Missing values in sample matrix",
        "",
        f"- Missing values: {missing_values:,}",
        f"- Total numeric cells: {total_values:,}",
        f"- Missing percentage: {missing_percent:.4f}%",
        "",
        "## Value range",
        "",
        f"- Minimum value: {minimum_value}",
        f"- Maximum value: {maximum_value}",
        "",
        "## First 10 sample identifiers",
        "",
        "```text",
        "\n".join(
            str(column)
            for column in sample_columns[:10]
        ),
        "```",
        "",
        "## First five rows",
        "",
        methylation.head().to_markdown(
            index=False
        ),
        "",
    ]

    return "\n".join(summary)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Run the methylation data preparation step."""

    print("Checking methylation files...")

    check_file(METHYLATION_450_FILE)
    print(
        f"[OK] {METHYLATION_450_FILE.name}"
    )

    check_file(METHYLATION_27_450_FILE)
    print(
        f"[OK] {METHYLATION_27_450_FILE.name}"
    )

    print(
        "\nLoading HM450 methylation dataset..."
    )

    methylation = load_methylation_data(
        METHYLATION_450_FILE
    )

    sample_columns = identify_sample_columns(
        methylation
    )

    metadata_columns = [
        column
        for column in methylation.columns
        if column not in sample_columns
    ]

    print(
        f"Rows/features: "
        f"{methylation.shape[0]:,}"
    )

    print(
        f"Total columns: "
        f"{methylation.shape[1]:,}"
    )

    print(
        f"Sample columns: "
        f"{len(sample_columns):,}"
    )

    print(
        f"Metadata columns: "
        f"{len(metadata_columns):,}"
    )

    print("\nMetadata column names:")
    print(metadata_columns)

    print("\nFirst 10 sample identifiers:")
    print(sample_columns[:10])

    print("\nCreating methylation report...")

    report = create_methylation_summary(
        methylation=methylation,
        file_path=METHYLATION_450_FILE,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")
    print(f"Report saved to:\n{REPORT_FILE}")


if __name__ == "__main__":
    main()