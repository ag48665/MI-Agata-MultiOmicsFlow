from pathlib import Path

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

INPUT_FILE = (
    PROCESSED_DIR
    / "clinical_genomic_epigenomic_merged.csv"
)

OUTPUT_FILE = (
    PROCESSED_DIR
    / "ml_epigenomic_feature_table.csv"
)

VARIANCE_TABLE_FILE = (
    TABLES_DIR
    / "methylation_feature_variance.csv"
)

SELECTED_FEATURES_FILE = (
    TABLES_DIR
    / "selected_methylation_features.csv"
)

REPORT_FILE = (
    REPORTS_DIR
    / "epigenomic_feature_selection.md"
)


# =====================================
# SETTINGS
# =====================================

TOP_N_METHYLATION_FEATURES = 500


# =====================================
# FUNCTIONS
# =====================================

def load_dataset(
    file_path: Path,
) -> pd.DataFrame:
    """Load the integrated multi-omics dataset."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
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


def identify_non_methylation_columns(
    df: pd.DataFrame,
    methylation_columns: list[str],
) -> list[str]:
    """Identify clinical, genomic and target columns."""

    return [
        column
        for column in df.columns
        if column not in methylation_columns
    ]


def calculate_feature_variance(
    df: pd.DataFrame,
    methylation_columns: list[str],
) -> pd.DataFrame:
    """Calculate variance for every methylation feature."""

    methylation_data = df[
        methylation_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    variance = methylation_data.var(
        axis=0,
        skipna=True,
    )

    variance_table = (
        variance
        .rename("variance")
        .reset_index()
        .rename(
            columns={
                "index": "feature",
            }
        )
        .sort_values(
            "variance",
            ascending=False,
        )
        .reset_index(drop=True)
    )

    variance_table["rank"] = (
        variance_table.index + 1
    )

    variance_table["missing_values"] = [
        int(
            methylation_data[
                feature
            ].isna().sum()
        )
        for feature in variance_table[
            "feature"
        ]
    ]

    variance_table["missing_percentage"] = [
        (
            methylation_data[
                feature
            ].isna().mean()
            * 100
        )
        for feature in variance_table[
            "feature"
        ]
    ]

    return variance_table[
        [
            "rank",
            "feature",
            "variance",
            "missing_values",
            "missing_percentage",
        ]
    ]


def select_top_features(
    variance_table: pd.DataFrame,
    top_n: int,
) -> list[str]:
    """Select the most variable methylation features."""

    valid_features = variance_table[
        variance_table["variance"].notna()
        & (variance_table["variance"] > 0)
    ].copy()

    selected_features = (
        valid_features
        .head(top_n)["feature"]
        .tolist()
    )

    if not selected_features:
        raise ValueError(
            "No valid methylation features "
            "were available for selection."
        )

    return selected_features


def build_feature_table(
    df: pd.DataFrame,
    non_methylation_columns: list[str],
    selected_methylation_features: list[str],
) -> pd.DataFrame:
    """Build the final AI-ready feature table."""

    selected_columns = (
        non_methylation_columns
        + selected_methylation_features
    )

    feature_table = df[
        selected_columns
    ].copy()

    if "patient_id" in feature_table.columns:
        feature_table["patient_id"] = (
            feature_table["patient_id"]
            .astype("string")
            .str.strip()
            .str.upper()
        )

    feature_table = (
        feature_table
        .drop_duplicates(
            subset=["patient_id"]
            if "patient_id" in feature_table.columns
            else None
        )
        .reset_index(drop=True)
    )

    return feature_table


def create_selected_features_table(
    variance_table: pd.DataFrame,
    selected_features: list[str],
) -> pd.DataFrame:
    """Create a table describing selected features."""

    selected_table = variance_table[
        variance_table["feature"].isin(
            selected_features
        )
    ].copy()

    selected_table = (
        selected_table
        .sort_values(
            "rank"
        )
        .reset_index(drop=True)
    )

    selected_table["selected"] = True

    return selected_table


def create_report(
    input_df: pd.DataFrame,
    feature_table: pd.DataFrame,
    all_methylation_columns: list[str],
    selected_features_table: pd.DataFrame,
    non_methylation_columns: list[str],
) -> str:
    """Create a Markdown feature-selection report."""

    top_20 = selected_features_table.head(20)

    top_feature_lines = [
        (
            f"{int(row.rank)}. `{row.feature}` — "
            f"variance: {row.variance:.6f}"
        )
        for row in top_20.itertuples(
            index=False
        )
    ]

    report_lines = [
        "# Epigenomic Feature Selection Report",
        "",
        "## Input dataset",
        "",
        f"- Patients: {input_df.shape[0]:,}",
        f"- Total input columns: {input_df.shape[1]:,}",
        (
            "- Total methylation features: "
            f"{len(all_methylation_columns):,}"
        ),
        (
            "- Non-methylation columns: "
            f"{len(non_methylation_columns):,}"
        ),
        "",
        "## Feature-selection method",
        "",
        (
            "- Methylation features were converted "
            "to numeric values."
        ),
        (
            "- Variance was calculated across patients "
            "for every methylation feature."
        ),
        (
            "- Features with missing or zero variance "
            "were excluded."
        ),
        (
            f"- The top "
            f"{len(selected_features_table):,} "
            "features with the highest variance "
            "were retained."
        ),
        "",
        "## Output feature table",
        "",
        f"- Patients: {feature_table.shape[0]:,}",
        f"- Final columns: {feature_table.shape[1]:,}",
        (
            "- Selected methylation features: "
            f"{len(selected_features_table):,}"
        ),
        "",
        "## Top 20 selected methylation features",
        "",
        *top_feature_lines,
        "",
        "## Output files",
        "",
        "- `data/processed/ml_epigenomic_feature_table.csv`",
        "- `outputs/tables/methylation_feature_variance.csv`",
        "- `outputs/tables/selected_methylation_features.csv`",
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Build the AI-ready epigenomic feature table."""

    print("Loading integrated dataset...")

    dataset = load_dataset(
        INPUT_FILE
    )

    print(
        f"Input rows: "
        f"{dataset.shape[0]:,}"
    )

    print(
        f"Input columns: "
        f"{dataset.shape[1]:,}"
    )

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

    non_methylation_columns = (
        identify_non_methylation_columns(
            dataset,
            methylation_columns,
        )
    )

    print(
        f"Methylation features detected: "
        f"{len(methylation_columns):,}"
    )

    print(
        f"Non-methylation columns detected: "
        f"{len(non_methylation_columns):,}"
    )

    print("\nCalculating feature variance...")

    variance_table = (
        calculate_feature_variance(
            dataset,
            methylation_columns,
        )
    )

    print(
        "\nSelecting top methylation features..."
    )

    selected_features = select_top_features(
        variance_table,
        TOP_N_METHYLATION_FEATURES,
    )

    print(
        f"Selected methylation features: "
        f"{len(selected_features):,}"
    )

    selected_features_table = (
        create_selected_features_table(
            variance_table,
            selected_features,
        )
    )

    print("\nBuilding AI-ready feature table...")

    feature_table = build_feature_table(
        dataset,
        non_methylation_columns,
        selected_features,
    )

    print(
        f"Final rows: "
        f"{feature_table.shape[0]:,}"
    )

    print(
        f"Final columns: "
        f"{feature_table.shape[1]:,}"
    )

    print("\nSaving outputs...")

    feature_table.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    variance_table.to_csv(
        VARIANCE_TABLE_FILE,
        index=False,
    )

    selected_features_table.to_csv(
        SELECTED_FEATURES_FILE,
        index=False,
    )

    report = create_report(
        input_df=dataset,
        feature_table=feature_table,
        all_methylation_columns=methylation_columns,
        selected_features_table=selected_features_table,
        non_methylation_columns=non_methylation_columns,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")

    print(
        "AI-ready feature table saved to:\n"
        f"{OUTPUT_FILE}"
    )

    print(
        "\nVariance table saved to:\n"
        f"{VARIANCE_TABLE_FILE}"
    )

    print(
        "\nSelected features saved to:\n"
        f"{SELECTED_FEATURES_FILE}"
    )

    print(
        "\nReport saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()