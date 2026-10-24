from pathlib import Path

import pandas as pd


# =====================================
# ŚCIEŻKI
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"

TABLES_DIR.mkdir(parents=True, exist_ok=True)

CLINICAL_PATIENT_FILE = (
    PROCESSED_DIR / "clinical_patient_cleaned.csv"
)

CLINICAL_SAMPLE_FILE = (
    PROCESSED_DIR / "clinical_sample_cleaned.csv"
)

GENOMICS_FILE = (
    PROCESSED_DIR / "genomics_cleaned.csv"
)

MERGED_OUTPUT_FILE = (
    PROCESSED_DIR / "clinical_genomics_merged.csv"
)

MATCHING_SUMMARY_FILE = (
    TABLES_DIR / "id_matching_summary.csv"
)


# =====================================
# FUNKCJE
# =====================================

def load_csv(file_path: Path) -> pd.DataFrame:
    """Wczytuje plik CSV i sprawdza, czy istnieje."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False
    )


def standardize_id_column(
    df: pd.DataFrame,
    column_name: str
) -> pd.DataFrame:
    """Czyści i standaryzuje kolumnę identyfikatora."""

    if column_name not in df.columns:
        raise KeyError(
            f"Nie znaleziono kolumny '{column_name}'."
        )

    df = df.copy()

    df[column_name] = (
        df[column_name]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return df


def create_mutation_features(
    genomics: pd.DataFrame
) -> pd.DataFrame:
    """
    Tworzy jeden rekord na pacjenta z podstawowymi
    cechami genomowymi.
    """

    mutation_features = (
        genomics.groupby("patient_id")
        .agg(
            mutation_count=("hugo_symbol", "size"),
            unique_mutated_genes=(
                "hugo_symbol",
                "nunique"
            ),
            unique_variant_classes=(
                "variant_classification",
                "nunique"
            ),
            genomic_sample_count=(
                "sample_id",
                "nunique"
            )
        )
        .reset_index()
    )

    return mutation_features


def create_matching_summary(
    clinical: pd.DataFrame,
    genomic_features: pd.DataFrame,
    merged: pd.DataFrame
) -> pd.DataFrame:
    """Tworzy raport dopasowania identyfikatorów."""

    clinical_ids = set(
        clinical["patient_id"].dropna()
    )

    genomic_ids = set(
        genomic_features["patient_id"].dropna()
    )

    matched_ids = clinical_ids.intersection(
        genomic_ids
    )

    unmatched_clinical = clinical_ids.difference(
        genomic_ids
    )

    unmatched_genomic = genomic_ids.difference(
        clinical_ids
    )

    match_percentage = (
        len(matched_ids) / len(clinical_ids) * 100
        if clinical_ids
        else 0
    )

    return pd.DataFrame(
        {
            "metric": [
                "clinical_records",
                "unique_clinical_patients",
                "unique_genomic_patients",
                "matched_patients",
                "unmatched_clinical_patients",
                "unmatched_genomic_patients",
                "merged_records",
                "match_percentage"
            ],
            "value": [
                len(clinical),
                len(clinical_ids),
                len(genomic_ids),
                len(matched_ids),
                len(unmatched_clinical),
                len(unmatched_genomic),
                len(merged),
                round(match_percentage, 2)
            ]
        }
    )


# =====================================
# MAIN
# =====================================

def main() -> None:
    print("Loading cleaned clinical data...")

    clinical_patient = load_csv(
        CLINICAL_PATIENT_FILE
    )

    clinical_sample = load_csv(
        CLINICAL_SAMPLE_FILE
    )

    print("Loading cleaned genomic data...")

    genomics = load_csv(
        GENOMICS_FILE
    )

    print("Standardizing identifiers...")

    clinical_patient = standardize_id_column(
        clinical_patient,
        "patient_id"
    )

    clinical_sample = standardize_id_column(
        clinical_sample,
        "patient_id"
    )

    genomics = standardize_id_column(
        genomics,
        "patient_id"
    )

    print("Preparing clinical table...")

    clinical = clinical_patient.merge(
        clinical_sample,
        on="patient_id",
        how="left",
        suffixes=("_patient", "_sample")
    )

    print("Creating patient-level genomic features...")

    genomic_features = create_mutation_features(
        genomics
    )

    print("Matching clinical and genomic data...")

    merged = clinical.merge(
        genomic_features,
        on="patient_id",
        how="inner"
    )

    matching_summary = create_matching_summary(
        clinical=clinical,
        genomic_features=genomic_features,
        merged=merged
    )

    print("Saving outputs...")

    merged.to_csv(
        MERGED_OUTPUT_FILE,
        index=False
    )

    matching_summary.to_csv(
        MATCHING_SUMMARY_FILE,
        index=False
    )

    print("\nFinished successfully.")
    print(
        f"Clinical rows: {len(clinical):,}"
    )
    print(
        f"Genomic patients: "
        f"{len(genomic_features):,}"
    )
    print(
        f"Matched rows: {len(merged):,}"
    )

    print("\nSaved files:")
    print(MERGED_OUTPUT_FILE)
    print(MATCHING_SUMMARY_FILE)


if __name__ == "__main__":
    main()