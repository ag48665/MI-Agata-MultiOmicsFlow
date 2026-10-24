from pathlib import Path

import pandas as pd


# =========================
# Ścieżki projektu
# =========================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)


PATIENT_FILE = RAW_DATA_DIR / "data_clinical_patient.txt"
SAMPLE_FILE = RAW_DATA_DIR / "data_clinical_sample.txt"


# =========================
# Funkcje
# =========================

def load_cbioportal_file(file_path: Path) -> pd.DataFrame:
    """Wczytuje plik TSV z cBioPortal."""
    return pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        low_memory=False
    )


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Ujednolica nazwy kolumn."""

    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(" ", "_")
        .str.replace("-", "_")
    )

    return df


def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Usuwa zduplikowane rekordy."""

    before = len(df)

    df = df.drop_duplicates()

    after = len(df)

    print(f"Usunięto {before-after} duplikatów.")

    return df


def create_missingness_table(df: pd.DataFrame, dataset_name: str):

    missing = df.isna().sum()

    summary = pd.DataFrame({

        "column": missing.index,
        "missing_values": missing.values,
        "missing_percent": (
            missing.values / len(df) * 100
        ).round(2),
        "dataset": dataset_name
    })

    return summary.sort_values(
        "missing_percent",
        ascending=False
    )


# =========================
# Main
# =========================

def main():

    print("Ładowanie danych...")

    patient = load_cbioportal_file(PATIENT_FILE)
    sample = load_cbioportal_file(SAMPLE_FILE)

    print("Standaryzacja nazw kolumn...")

    patient = standardize_columns(patient)
    sample = standardize_columns(sample)

    print("Usuwanie duplikatów...")

    patient = remove_duplicates(patient)
    sample = remove_duplicates(sample)

    print("Tworzenie raportu brakujących danych...")

    patient_missing = create_missingness_table(
        patient,
        "clinical_patient"
    )

    sample_missing = create_missingness_table(
        sample,
        "clinical_sample"
    )

    missing_summary = pd.concat(
        [patient_missing, sample_missing],
        ignore_index=True
    )

    print("Zapisywanie plików...")

    patient.to_csv(
        PROCESSED_DIR / "clinical_patient_cleaned.csv",
        index=False
    )

    sample.to_csv(
        PROCESSED_DIR / "clinical_sample_cleaned.csv",
        index=False
    )

    missing_summary.to_csv(
        TABLES_DIR / "clinical_missingness_summary.csv",
        index=False
    )

    print("\nGotowe!")

    print(
        "\nPliki zapisano:"
    )

    print(
        PROCESSED_DIR /
        "clinical_patient_cleaned.csv"
    )

    print(
        PROCESSED_DIR /
        "clinical_sample_cleaned.csv"
    )

    print(
        TABLES_DIR /
        "clinical_missingness_summary.csv"
    )


if __name__ == "__main__":
    main()