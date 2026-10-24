from pathlib import Path

import pandas as pd


# Główny folder projektu
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Folder z surowymi danymi
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA"

# Pliki wejściowe
CLINICAL_PATIENT_FILE = RAW_DATA_DIR / "data_clinical_patient.txt"
CLINICAL_SAMPLE_FILE = RAW_DATA_DIR / "data_clinical_sample.txt"
MUTATION_FILE = RAW_DATA_DIR / "data_mutations.txt"

# Folder na raporty
REPORTS_DIR = PROJECT_ROOT / "outputs" / "reports"
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

REPORT_FILE = REPORTS_DIR / "data_overview.md"


def check_file(file_path: Path) -> None:
    """Sprawdza, czy wymagany plik istnieje."""
    if not file_path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{file_path}\n"
            "Sprawdź, czy plik został skopiowany do data/raw/TCGA-BRCA."
        )


def load_cbioportal_file(file_path: Path) -> pd.DataFrame:
    """
    Wczytuje plik cBioPortal zapisany jako tabela TSV.

    Wiersze zaczynające się od # zawierają metadane,
    dlatego są pomijane podczas wczytywania.
    """
    return pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        low_memory=False,
    )


def create_dataset_summary(name: str, dataframe: pd.DataFrame) -> str:
    """Tworzy tekstowe podsumowanie zbioru danych."""
    missing_values = dataframe.isna().sum()
    missing_values = missing_values[missing_values > 0]
    missing_values = missing_values.sort_values(ascending=False)

    summary = [
        f"## {name}",
        "",
        f"- Liczba wierszy: {dataframe.shape[0]}",
        f"- Liczba kolumn: {dataframe.shape[1]}",
        f"- Liczba duplikatów: {dataframe.duplicated().sum()}",
        "",
        "### Nazwy kolumn",
        "",
        "```text",
        "\n".join(dataframe.columns.astype(str)),
        "```",
        "",
        "### Brakujące wartości",
        "",
    ]

    if missing_values.empty:
        summary.append("Nie znaleziono brakujących wartości.")
    else:
        summary.extend(
            [
                "```text",
                missing_values.to_string(),
                "```",
            ]
        )

    summary.extend(
        [
            "",
            "### Pierwsze 5 wierszy",
            "",
            dataframe.head().to_markdown(index=False),
            "",
        ]
    )

    return "\n".join(summary)


def main() -> None:
    required_files = [
        CLINICAL_PATIENT_FILE,
        CLINICAL_SAMPLE_FILE,
        MUTATION_FILE,
    ]

    print("Sprawdzanie wymaganych plików...")

    for file_path in required_files:
        check_file(file_path)
        print(f"[OK] {file_path.name}")

    print("\nWczytywanie danych...")

    clinical_patient = load_cbioportal_file(CLINICAL_PATIENT_FILE)
    clinical_sample = load_cbioportal_file(CLINICAL_SAMPLE_FILE)
    mutations = load_cbioportal_file(MUTATION_FILE)

    datasets = {
        "Clinical patient data": clinical_patient,
        "Clinical sample data": clinical_sample,
        "Mutation data": mutations,
    }

    report_sections = [
        "# TCGA-BRCA Data Overview",
        "",
        "Źródło danych: cBioPortal, TCGA PanCancer Atlas.",
        "",
    ]

    for name, dataframe in datasets.items():
        print(f"\n{name}")
        print(f"Wymiary: {dataframe.shape[0]} wierszy × {dataframe.shape[1]} kolumn")
        print(dataframe.head())

        report_sections.append(
            create_dataset_summary(name, dataframe)
        )

    REPORT_FILE.write_text(
        "\n".join(report_sections),
        encoding="utf-8",
    )

    print("\nZakończono przygotowanie danych.")
    print(f"Raport zapisano w:\n{REPORT_FILE}")


if __name__ == "__main__":
    main()