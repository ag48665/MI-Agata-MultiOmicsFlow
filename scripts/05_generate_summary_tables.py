from pathlib import Path

import pandas as pd


# =====================================
# ŚCIEŻKI
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"

TABLES_DIR.mkdir(parents=True, exist_ok=True)

MERGED_FILE = PROCESSED_DIR / "clinical_genomics_merged.csv"
GENOMICS_FILE = PROCESSED_DIR / "genomics_cleaned.csv"

CLINICAL_SUMMARY_FILE = (
    TABLES_DIR / "clinical_genomic_summary.csv"
)

MUTATION_BURDEN_FILE = (
    TABLES_DIR / "mutation_burden_summary.csv"
)

SUBTYPE_SUMMARY_FILE = (
    TABLES_DIR / "cancer_subtype_summary.csv"
)

SURVIVAL_SUMMARY_FILE = (
    TABLES_DIR / "survival_status_summary.csv"
)

TOP_MUTATED_GENES_FILE = (
    TABLES_DIR / "top_mutated_genes.csv"
)


# =====================================
# FUNKCJE
# =====================================

def load_csv(file_path: Path) -> pd.DataFrame:
    """Wczytuje plik CSV i sprawdza jego istnienie."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False
    )


def find_first_existing_column(
    df: pd.DataFrame,
    candidates: list[str]
) -> str | None:
    """Zwraca pierwszą istniejącą kolumnę z listy."""

    for column in candidates:
        if column in df.columns:
            return column

    return None


def create_numeric_summary(
    df: pd.DataFrame,
    columns: list[str]
) -> pd.DataFrame:
    """Tworzy statystyki opisowe dla kolumn numerycznych."""

    available_columns = [
        column
        for column in columns
        if column in df.columns
    ]

    if not available_columns:
        return pd.DataFrame(
            columns=[
                "variable",
                "count",
                "mean",
                "std",
                "min",
                "median",
                "max"
            ]
        )

    summaries = []

    for column in available_columns:
        values = pd.to_numeric(
            df[column],
            errors="coerce"
        )

        summaries.append(
            {
                "variable": column,
                "count": int(values.count()),
                "mean": round(values.mean(), 2),
                "std": round(values.std(), 2),
                "min": round(values.min(), 2),
                "median": round(values.median(), 2),
                "max": round(values.max(), 2)
            }
        )

    return pd.DataFrame(summaries)


def create_categorical_summary(
    df: pd.DataFrame,
    column_name: str
) -> pd.DataFrame:
    """Tworzy podsumowanie wartości kategorii."""

    if column_name not in df.columns:
        return pd.DataFrame(
            columns=[
                column_name,
                "count",
                "percent"
            ]
        )

    summary = (
        df[column_name]
        .fillna("Missing")
        .astype(str)
        .value_counts(dropna=False)
        .reset_index()
    )

    summary.columns = [
        column_name,
        "count"
    ]

    summary["percent"] = (
        summary["count"]
        / len(df)
        * 100
    ).round(2)

    return summary


def create_top_mutated_genes(
    genomics: pd.DataFrame,
    top_n: int = 20
) -> pd.DataFrame:
    """Tworzy ranking najczęściej mutowanych genów."""

    required_columns = {
        "hugo_symbol",
        "patient_id",
        "sample_id"
    }

    missing = required_columns.difference(
        genomics.columns
    )

    if missing:
        raise KeyError(
            f"Brak wymaganych kolumn: {sorted(missing)}"
        )

    return (
        genomics.groupby("hugo_symbol")
        .agg(
            mutation_count=("hugo_symbol", "size"),
            patient_count=("patient_id", "nunique"),
            sample_count=("sample_id", "nunique")
        )
        .reset_index()
        .rename(columns={"hugo_symbol": "gene"})
        .sort_values(
            ["patient_count", "mutation_count"],
            ascending=False
        )
        .head(top_n)
        .reset_index(drop=True)
    )


# =====================================
# MAIN
# =====================================

def main() -> None:
    print("Loading merged clinical-genomic data...")

    merged = load_csv(MERGED_FILE)

    print("Loading cleaned genomic data...")

    genomics = load_csv(GENOMICS_FILE)

    print(
        f"Merged dataset: "
        f"{merged.shape[0]:,} rows x "
        f"{merged.shape[1]:,} columns"
    )

    print("Identifying clinical columns...")

    age_column = find_first_existing_column(
        merged,
        [
            "age",
            "age_at_diagnosis",
            "age_patient"
        ]
    )

    subtype_column = find_first_existing_column(
        merged,
        [
            "subtype",
            "subtype_patient",
            "cancer_type_detailed",
            "cancer_type_detailed_sample"
        ]
    )

    survival_column = find_first_existing_column(
        merged,
        [
            "os_status",
            "vital_status",
            "survival_status",
            "os_status_patient"
        ]
    )

    print("Generating clinical-genomic summary...")

    numeric_candidates = [
        column
        for column in [
            age_column,
            "mutation_count",
            "unique_mutated_genes",
            "unique_variant_classes",
            "genomic_sample_count",
            "os_months",
            "dfs_months",
            "pfs_months",
            "dss_months"
        ]
        if column is not None
    ]

    clinical_summary = create_numeric_summary(
        merged,
        numeric_candidates
    )

    print("Generating mutation burden summary...")

    mutation_burden_summary = create_numeric_summary(
        merged,
        [
            "mutation_count",
            "unique_mutated_genes",
            "unique_variant_classes"
        ]
    )

    print("Generating subtype summary...")

    if subtype_column is not None:
        subtype_summary = create_categorical_summary(
            merged,
            subtype_column
        )
    else:
        subtype_summary = pd.DataFrame(
            columns=[
                "subtype",
                "count",
                "percent"
            ]
        )

    print("Generating survival summary...")

    if survival_column is not None:
        survival_summary = create_categorical_summary(
            merged,
            survival_column
        )
    else:
        survival_summary = pd.DataFrame(
            columns=[
                "survival_status",
                "count",
                "percent"
            ]
        )

    print("Generating top mutated genes...")

    top_genes = create_top_mutated_genes(
        genomics,
        top_n=20
    )

    print("Saving outputs...")

    clinical_summary.to_csv(
        CLINICAL_SUMMARY_FILE,
        index=False
    )

    mutation_burden_summary.to_csv(
        MUTATION_BURDEN_FILE,
        index=False
    )

    subtype_summary.to_csv(
        SUBTYPE_SUMMARY_FILE,
        index=False
    )

    survival_summary.to_csv(
        SURVIVAL_SUMMARY_FILE,
        index=False
    )

    top_genes.to_csv(
        TOP_MUTATED_GENES_FILE,
        index=False
    )

    print("\nFinished successfully.")

    print("\nDetected columns:")
    print(f"Age column: {age_column}")
    print(f"Subtype column: {subtype_column}")
    print(f"Survival column: {survival_column}")

    print("\nSaved files:")
    print(CLINICAL_SUMMARY_FILE)
    print(MUTATION_BURDEN_FILE)
    print(SUBTYPE_SUMMARY_FILE)
    print(SURVIVAL_SUMMARY_FILE)
    print(TOP_MUTATED_GENES_FILE)


if __name__ == "__main__":
    main()