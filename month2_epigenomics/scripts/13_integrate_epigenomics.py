from pathlib import Path
import re

import pandas as pd


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
)

METHYLATION_DIR = (
    PROCESSED_DIR
    / "methylation"
)

TABLES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "tables"
)

REPORTS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
)

TABLES_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

REPORTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

METHYLATION_FILE = (
    METHYLATION_DIR
    / "hm450_clean.csv"
)

MAPPING_FILE = (
    METHYLATION_DIR
    / "methylation_patient_mapping.csv"
)

CLINICAL_GENOMICS_FILE = (
    PROCESSED_DIR
    / "clinical_genomics_merged.csv"
)

METHYLATION_PATIENT_FILE = (
    METHYLATION_DIR
    / "methylation_patient_matrix.csv"
)

INTEGRATED_OUTPUT_FILE = (
    PROCESSED_DIR
    / "clinical_genomic_epigenomic_merged.csv"
)

FEATURE_ANNOTATION_FILE = (
    TABLES_DIR
    / "methylation_feature_annotation.csv"
)

INTEGRATION_SUMMARY_FILE = (
    TABLES_DIR
    / "epigenomic_integration_summary.csv"
)

REPORT_FILE = (
    REPORTS_DIR
    / "epigenomic_integration_report.md"
)


# =====================================
# FUNCTIONS
# =====================================

def load_csv(file_path: Path) -> pd.DataFrame:
    """Load a CSV file and verify that it exists."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def standardize_identifier(
    values: pd.Series,
) -> pd.Series:
    """Standardize TCGA identifiers."""

    return (
        values
        .astype("string")
        .str.strip()
        .str.upper()
    )


def identify_sample_columns(
    methylation: pd.DataFrame,
) -> list[str]:
    """Identify TCGA sample columns."""

    return [
        column
        for column in methylation.columns
        if str(column).upper().startswith("TCGA-")
    ]


def identify_metadata_columns(
    methylation: pd.DataFrame,
    sample_columns: list[str],
) -> list[str]:
    """Return columns containing feature annotations."""

    return [
        column
        for column in methylation.columns
        if column not in sample_columns
    ]


def find_feature_identifier_column(
    methylation: pd.DataFrame,
) -> str:
    """Find the best column for naming methylation features."""

    candidates = [
        "ENTIITY_STABLE_ID",
        "ENTITY_STABLE_ID",
        "NAME",
        "TRANSCRIPT_ID",
        "DESCRIPTION",
    ]

    for column in candidates:
        if column in methylation.columns:
            return column

    raise KeyError(
        "No suitable methylation feature identifier "
        "column was detected."
    )


def sanitize_feature_name(value: object) -> str:
    """Convert a feature identifier into a safe column name."""

    text = str(value).strip()

    if not text or text.upper() in {"NAN", "NONE", "<NA>"}:
        text = "unknown"

    text = re.sub(
        r"[^A-Za-z0-9_]+",
        "_",
        text,
    )

    text = text.strip("_")

    return text or "unknown"


def make_unique_feature_names(
    methylation: pd.DataFrame,
    identifier_column: str,
) -> list[str]:
    """Create unique feature names for the wide patient matrix."""

    counts: dict[str, int] = {}
    feature_names: list[str] = []

    for row_number, value in enumerate(
        methylation[identifier_column],
        start=1,
    ):
        base_name = sanitize_feature_name(value)

        if base_name == "unknown":
            base_name = f"feature_{row_number}"

        base_name = f"meth_{base_name}"

        count = counts.get(base_name, 0)

        if count == 0:
            unique_name = base_name
        else:
            unique_name = f"{base_name}_{count + 1}"

        counts[base_name] = count + 1
        feature_names.append(unique_name)

    return feature_names


def create_feature_annotation(
    methylation: pd.DataFrame,
    metadata_columns: list[str],
    feature_names: list[str],
) -> pd.DataFrame:
    """Create a table linking new feature names to annotations."""

    annotation = methylation[
        metadata_columns
    ].copy()

    annotation.insert(
        0,
        "methylation_feature",
        feature_names,
    )

    return annotation


def create_patient_methylation_matrix(
    methylation: pd.DataFrame,
    mapping: pd.DataFrame,
    sample_columns: list[str],
    feature_names: list[str],
) -> pd.DataFrame:
    """
    Convert methylation data from:

    features x samples

    into:

    patients x features
    """

    sample_matrix = methylation[
        sample_columns
    ].copy()

    sample_matrix = sample_matrix.apply(
        pd.to_numeric,
        errors="coerce",
    )

    sample_matrix.index = feature_names

    patient_matrix = (
        sample_matrix
        .transpose()
        .reset_index()
        .rename(
            columns={
                "index": "sample_id",
            }
        )
    )

    patient_matrix["sample_id"] = (
        standardize_identifier(
            patient_matrix["sample_id"]
        )
    )

    mapping = mapping.copy()

    mapping["sample_id"] = standardize_identifier(
        mapping["sample_id"]
    )

    mapping["patient_id"] = standardize_identifier(
        mapping["patient_id"]
    )

    matched_matrix = patient_matrix.merge(
        mapping[
            [
                "sample_id",
                "patient_id",
            ]
        ],
        on="sample_id",
        how="inner",
        validate="one_to_one",
    )

    methylation_features = [
        column
        for column in matched_matrix.columns
        if column.startswith("meth_")
    ]

    patient_level_matrix = (
        matched_matrix
        .groupby(
            "patient_id",
            as_index=False,
        )[methylation_features]
        .mean()
    )

    return patient_level_matrix


def integrate_datasets(
    clinical_genomics: pd.DataFrame,
    methylation_patients: pd.DataFrame,
) -> tuple[pd.DataFrame, int]:
    """Integrate clinical-genomic and methylation datasets."""

    if "patient_id" not in clinical_genomics.columns:
        raise KeyError(
            "Column 'patient_id' was not found in "
            "clinical_genomics_merged.csv."
        )

    clinical_genomics = clinical_genomics.copy()

    clinical_genomics["patient_id"] = (
        standardize_identifier(
            clinical_genomics["patient_id"]
        )
    )

    duplicate_clinical_rows = int(
        clinical_genomics.duplicated(
            subset=["patient_id"]
        ).sum()
    )

    clinical_genomics = (
        clinical_genomics
        .drop_duplicates(
            subset=["patient_id"],
            keep="first",
        )
        .copy()
    )

    integrated = clinical_genomics.merge(
        methylation_patients,
        on="patient_id",
        how="inner",
        validate="one_to_one",
    )

    return integrated, duplicate_clinical_rows


def create_integration_summary(
    methylation: pd.DataFrame,
    mapping: pd.DataFrame,
    clinical_genomics: pd.DataFrame,
    methylation_patients: pd.DataFrame,
    integrated: pd.DataFrame,
    duplicate_clinical_rows: int,
    feature_count: int,
) -> pd.DataFrame:
    """Create integration summary statistics."""

    clinical_patient_count = (
        clinical_genomics["patient_id"]
        .astype("string")
        .str.strip()
        .str.upper()
        .nunique()
    )

    return pd.DataFrame(
        {
            "metric": [
                "input_methylation_features",
                "input_methylation_samples",
                "mapped_methylation_samples",
                "unique_methylation_patients",
                "unique_clinical_genomic_patients",
                "matched_integrated_patients",
                "methylation_features_added",
                "duplicate_clinical_patient_rows_removed",
                "integrated_total_columns",
            ],
            "value": [
                methylation.shape[0],
                len(
                    identify_sample_columns(
                        methylation
                    )
                ),
                len(mapping),
                methylation_patients[
                    "patient_id"
                ].nunique(),
                clinical_patient_count,
                len(integrated),
                feature_count,
                duplicate_clinical_rows,
                integrated.shape[1],
            ],
        }
    )


def create_report(
    summary: pd.DataFrame,
    integrated: pd.DataFrame,
) -> str:
    """Create a Markdown integration report."""

    summary_lines = [
        (
            f"- {row.metric}: "
            f"{row.value:,}"
        )
        for row in summary.itertuples(
            index=False
        )
    ]

    report_lines = [
        "# Clinical–Genomic–Epigenomic Integration Report",
        "",
        "## Integration summary",
        "",
        *summary_lines,
        "",
        "## Integrated dataset",
        "",
        (
            "- Output file: "
            "`clinical_genomic_epigenomic_merged.csv`"
        ),
        (
            f"- Rows/patients: "
            f"{integrated.shape[0]:,}"
        ),
        (
            f"- Total columns: "
            f"{integrated.shape[1]:,}"
        ),
        "",
        "## Processing steps",
        "",
        (
            "- Converted the methylation matrix from "
            "features × samples to patients × features."
        ),
        (
            "- Converted TCGA sample identifiers to "
            "patient identifiers using the Month 2 mapping."
        ),
        (
            "- Averaged methylation values when more than "
            "one sample represented the same patient."
        ),
        (
            "- Integrated methylation features with the "
            "Month 1 clinical-genomic dataset."
        ),
        (
            "- Retained only patients available in both "
            "clinical-genomic and methylation datasets."
        ),
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Integrate clinical, genomic and epigenomic data."""

    print("Loading methylation data...")

    methylation = load_csv(
        METHYLATION_FILE
    )

    print("Loading methylation patient mapping...")

    mapping = load_csv(
        MAPPING_FILE
    )

    print("Loading clinical-genomic data...")

    clinical_genomics = load_csv(
        CLINICAL_GENOMICS_FILE
    )

    sample_columns = identify_sample_columns(
        methylation
    )

    if not sample_columns:
        raise ValueError(
            "No TCGA methylation sample columns "
            "were detected."
        )

    metadata_columns = identify_metadata_columns(
        methylation,
        sample_columns,
    )

    identifier_column = (
        find_feature_identifier_column(
            methylation
        )
    )

    print(
        f"Methylation features: "
        f"{methylation.shape[0]:,}"
    )

    print(
        f"Methylation samples: "
        f"{len(sample_columns):,}"
    )

    print(
        f"Feature identifier column: "
        f"{identifier_column}"
    )

    print("\nCreating unique methylation features...")

    feature_names = make_unique_feature_names(
        methylation,
        identifier_column,
    )

    feature_annotation = create_feature_annotation(
        methylation,
        metadata_columns,
        feature_names,
    )

    print(
        "\nCreating patient-level methylation matrix..."
    )

    methylation_patients = (
        create_patient_methylation_matrix(
            methylation=methylation,
            mapping=mapping,
            sample_columns=sample_columns,
            feature_names=feature_names,
        )
    )

    print(
        f"Patient-level methylation rows: "
        f"{methylation_patients.shape[0]:,}"
    )

    print(
        f"Methylation feature columns: "
        f"{methylation_patients.shape[1] - 1:,}"
    )

    print(
        "\nIntegrating clinical, genomic "
        "and methylation data..."
    )

    integrated, duplicate_clinical_rows = (
        integrate_datasets(
            clinical_genomics,
            methylation_patients,
        )
    )

    print(
        f"Integrated patients: "
        f"{integrated.shape[0]:,}"
    )

    print(
        f"Integrated columns: "
        f"{integrated.shape[1]:,}"
    )

    print("\nCreating integration summary...")

    summary = create_integration_summary(
        methylation=methylation,
        mapping=mapping,
        clinical_genomics=clinical_genomics,
        methylation_patients=methylation_patients,
        integrated=integrated,
        duplicate_clinical_rows=duplicate_clinical_rows,
        feature_count=len(feature_names),
    )

    print("\nSaving outputs...")

    methylation_patients.to_csv(
        METHYLATION_PATIENT_FILE,
        index=False,
    )

    integrated.to_csv(
        INTEGRATED_OUTPUT_FILE,
        index=False,
    )

    feature_annotation.to_csv(
        FEATURE_ANNOTATION_FILE,
        index=False,
    )

    summary.to_csv(
        INTEGRATION_SUMMARY_FILE,
        index=False,
    )

    report = create_report(
        summary,
        integrated,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")

    print(
        "Patient methylation matrix saved to:\n"
        f"{METHYLATION_PATIENT_FILE}"
    )

    print(
        "\nIntegrated dataset saved to:\n"
        f"{INTEGRATED_OUTPUT_FILE}"
    )

    print(
        "\nFeature annotation saved to:\n"
        f"{FEATURE_ANNOTATION_FILE}"
    )

    print(
        "\nIntegration summary saved to:\n"
        f"{INTEGRATION_SUMMARY_FILE}"
    )

    print(
        "\nReport saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()