
from pathlib import Path

import pandas as pd


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "TCGA-BRCA"
)

PROCESSED_DATA_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "methylation"
)

REPORTS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
)

PROCESSED_DATA_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

REPORTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

INPUT_FILE = (
    RAW_DATA_DIR
    / "data_methylation_hm450.txt"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR
    / "hm450_clean.csv"
)

QC_REPORT = (
    REPORTS_DIR
    / "methylation_qc_report.md"
)


# =====================================
# FUNCTIONS
# =====================================

def load_data(file_path: Path) -> pd.DataFrame:
    """Load the methylation dataset."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        low_memory=False,
    )


def identify_sample_columns(
    df: pd.DataFrame,
) -> list[str]:
    """Identify TCGA sample columns."""

    return [
        column
        for column in df.columns
        if str(column).upper().startswith("TCGA-")
    ]


def quality_control(
    methylation: pd.DataFrame,
    sample_columns: list[str],
) -> dict:
    """Calculate basic methylation QC statistics."""

    sample_data = methylation[
        sample_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    total_values = int(
        sample_data.shape[0]
        * sample_data.shape[1]
    )

    missing_values = int(
        sample_data.isna().sum().sum()
    )

    missing_percentage = (
        missing_values / total_values * 100
        if total_values
        else 0
    )

    duplicated_rows = int(
        methylation.duplicated().sum()
    )

    values_below_zero = int(
        (sample_data < 0).sum().sum()
    )

    values_above_one = int(
        (sample_data > 1).sum().sum()
    )

    minimum_value = sample_data.min().min()
    maximum_value = sample_data.max().max()

    return {
        "total_values": total_values,
        "missing_values": missing_values,
        "missing_percentage": missing_percentage,
        "duplicated_rows": duplicated_rows,
        "values_below_zero": values_below_zero,
        "values_above_one": values_above_one,
        "minimum_value": minimum_value,
        "maximum_value": maximum_value,
    }


def clean_methylation_data(
    methylation: pd.DataFrame,
    sample_columns: list[str],
) -> pd.DataFrame:
    """Clean methylation values and remove duplicate rows."""

    cleaned = methylation.copy()

    cleaned[sample_columns] = cleaned[
        sample_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    cleaned = cleaned.drop_duplicates().copy()

    cleaned[sample_columns] = cleaned[
        sample_columns
    ].clip(
        lower=0,
        upper=1,
    )

    return cleaned


def create_qc_report(
    methylation: pd.DataFrame,
    cleaned: pd.DataFrame,
    sample_columns: list[str],
    metadata_columns: list[str],
    qc: dict,
) -> str:
    """Create a Markdown quality-control report."""

    report_lines = [
        "# TCGA-BRCA Methylation QC Report",
        "",
        "## Dataset dimensions",
        "",
        f"- Input rows/features: {methylation.shape[0]:,}",
        f"- Cleaned rows/features: {cleaned.shape[0]:,}",
        f"- Total columns: {methylation.shape[1]:,}",
        f"- Sample columns: {len(sample_columns):,}",
        f"- Metadata columns: {len(metadata_columns):,}",
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
        "## Quality-control results",
        "",
        f"- Total methylation values: {qc['total_values']:,}",
        f"- Missing values: {qc['missing_values']:,}",
        (
            "- Missing percentage: "
            f"{qc['missing_percentage']:.4f}%"
        ),
        f"- Duplicated rows: {qc['duplicated_rows']:,}",
        f"- Values below 0: {qc['values_below_zero']:,}",
        f"- Values above 1: {qc['values_above_one']:,}",
        f"- Minimum value: {qc['minimum_value']}",
        f"- Maximum value: {qc['maximum_value']}",
        "",
        "## Cleaning actions",
        "",
        "- Converted sample columns to numeric values.",
        "- Converted invalid entries to missing values.",
        "- Removed fully duplicated rows.",
        "- Restricted methylation values to the 0–1 range.",
        "- Retained missing values for later imputation.",
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Run methylation quality control and cleaning."""

    print("Loading methylation dataset...")

    methylation = load_data(
        INPUT_FILE
    )

    sample_columns = identify_sample_columns(
        methylation
    )

    if not sample_columns:
        raise ValueError(
            "No TCGA sample columns were detected."
        )

    metadata_columns = [
        column
        for column in methylation.columns
        if column not in sample_columns
    ]

    print(f"Rows: {methylation.shape[0]:,}")
    print(f"Columns: {methylation.shape[1]:,}")
    print(f"Samples: {len(sample_columns):,}")
    print(f"Metadata: {len(metadata_columns):,}")

    qc = quality_control(
        methylation,
        sample_columns,
    )

    print("\nQuality control")
    print("-" * 40)
    print(f"Missing values: {qc['missing_values']:,}")
    print(
        "Missing percentage: "
        f"{qc['missing_percentage']:.4f}%"
    )
    print(
        f"Duplicated rows: "
        f"{qc['duplicated_rows']:,}"
    )
    print(
        f"Values below 0: "
        f"{qc['values_below_zero']:,}"
    )
    print(
        f"Values above 1: "
        f"{qc['values_above_one']:,}"
    )
    print(
        f"Minimum value: "
        f"{qc['minimum_value']}"
    )
    print(
        f"Maximum value: "
        f"{qc['maximum_value']}"
    )

    print("\nCleaning methylation dataset...")

    cleaned = clean_methylation_data(
        methylation,
        sample_columns,
    )

    print("\nSaving cleaned dataset...")

    cleaned.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    print("\nCreating QC report...")

    report = create_qc_report(
        methylation=methylation,
        cleaned=cleaned,
        sample_columns=sample_columns,
        metadata_columns=metadata_columns,
        qc=qc,
    )

    QC_REPORT.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")
    print(
        f"Clean dataset saved to:\n{OUTPUT_FILE}"
    )
    print(
        f"QC report saved to:\n{QC_REPORT}"
    )


if __name__ == "__main__":
    main()
