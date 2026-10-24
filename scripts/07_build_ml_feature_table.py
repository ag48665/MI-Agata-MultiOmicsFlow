from pathlib import Path

import pandas as pd


# =====================================
# ŚCIEŻKI
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"

MERGED_FILE = (
    PROCESSED_DIR / "clinical_genomics_merged.csv"
)

GENOMICS_FILE = (
    PROCESSED_DIR / "genomics_cleaned.csv"
)

TOP_GENES_FILE = (
    TABLES_DIR / "top_mutated_genes.csv"
)

ML_OUTPUT_FILE = (
    PROCESSED_DIR / "ml_feature_table.csv"
)

ML_SUMMARY_FILE = (
    TABLES_DIR / "ml_feature_table_summary.csv"
)


# =====================================
# USTAWIENIA
# =====================================

TOP_N_GENES = 10

PREFERRED_GENES = [
    "PIK3CA",
    "TP53",
    "TTN",
    "CDH1",
    "GATA3",
    "MUC16",
    "KMT2C",
    "MAP3K1",
    "RYR2",
    "HMCN1"
]


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


def standardize_patient_id(
    df: pd.DataFrame
) -> pd.DataFrame:
    """Standaryzuje identyfikatory pacjentów."""

    if "patient_id" not in df.columns:
        raise KeyError(
            "Nie znaleziono kolumny 'patient_id'."
        )

    df = df.copy()

    df["patient_id"] = (
        df["patient_id"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return df


def find_first_existing_column(
    df: pd.DataFrame,
    candidates: list[str]
) -> str | None:
    """Zwraca pierwszą istniejącą kolumnę."""

    for column in candidates:
        if column in df.columns:
            return column

    return None


def choose_top_genes(
    top_genes: pd.DataFrame,
    genomics: pd.DataFrame,
    top_n: int
) -> list[str]:
    """Wybiera geny do kodowania jako cechy binarne."""

    available_genes = set(
        genomics["hugo_symbol"]
        .dropna()
        .astype(str)
        .str.upper()
    )

    selected_genes = [
        gene
        for gene in PREFERRED_GENES
        if gene in available_genes
    ]

    if "gene" in top_genes.columns:
        ranked_genes = (
            top_genes["gene"]
            .dropna()
            .astype(str)
            .str.upper()
            .tolist()
        )

        for gene in ranked_genes:
            if gene not in selected_genes:
                selected_genes.append(gene)

            if len(selected_genes) >= top_n:
                break

    return selected_genes[:top_n]


def create_gene_mutation_matrix(
    genomics: pd.DataFrame,
    selected_genes: list[str]
) -> pd.DataFrame:
    """
    Tworzy macierz binarną:
    1 = gen zmutowany u pacjenta,
    0 = brak mutacji.
    """

    subset = genomics[
        genomics["hugo_symbol"].isin(selected_genes)
    ][
        ["patient_id", "hugo_symbol"]
    ].drop_duplicates()

    subset["mutated"] = 1

    gene_matrix = subset.pivot_table(
        index="patient_id",
        columns="hugo_symbol",
        values="mutated",
        aggfunc="max",
        fill_value=0
    )

    gene_matrix = gene_matrix.reindex(
        columns=selected_genes,
        fill_value=0
    )

    gene_matrix.columns = [
        f"{gene}_mutated"
        for gene in gene_matrix.columns
    ]

    return gene_matrix.reset_index()


def create_outcome(
    df: pd.DataFrame,
    survival_column: str
) -> pd.Series:
    """
    Tworzy binarną zmienną docelową.

    0 = żyje / censored
    1 = zgon / deceased
    """

    values = (
        df[survival_column]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    outcome = pd.Series(
        pd.NA,
        index=df.index,
        dtype="Int64"
    )

    deceased_patterns = [
        "DECEASED",
        "DEAD",
        "1:DECEASED"
    ]

    living_patterns = [
        "LIVING",
        "ALIVE",
        "0:LIVING",
        "0:CENSORED"
    ]

    deceased_mask = values.apply(
        lambda value: any(
            pattern in str(value)
            for pattern in deceased_patterns
        )
    )

    living_mask = values.apply(
        lambda value: any(
            pattern in str(value)
            for pattern in living_patterns
        )
    )

    outcome.loc[deceased_mask] = 1
    outcome.loc[living_mask] = 0

    return outcome


def build_feature_table(
    merged: pd.DataFrame,
    gene_matrix: pd.DataFrame,
    selected_genes: list[str]
) -> pd.DataFrame:
    """Buduje końcową tabelę cech ML."""

    age_column = find_first_existing_column(
        merged,
        [
            "age",
            "age_at_diagnosis",
            "age_patient"
        ]
    )

    sex_column = find_first_existing_column(
        merged,
        [
            "sex",
            "gender",
            "sex_patient"
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

    if survival_column is None:
        raise KeyError(
            "Nie znaleziono kolumny przeżycia "
            "do utworzenia zmiennej outcome."
        )

    feature_table = pd.DataFrame(
        {
            "patient_id": merged["patient_id"]
        }
    )

    if age_column is not None:
        feature_table["age"] = pd.to_numeric(
            merged[age_column],
            errors="coerce"
        )

    if sex_column is not None:
        feature_table["sex"] = (
            merged[sex_column]
            .astype("string")
            .str.strip()
        )

    if subtype_column is not None:
        feature_table["subtype"] = (
            merged[subtype_column]
            .astype("string")
            .str.strip()
        )

    numeric_features = [
        "mutation_count",
        "unique_mutated_genes",
        "unique_variant_classes",
        "genomic_sample_count"
    ]

    for column in numeric_features:
        if column in merged.columns:
            feature_table[column] = pd.to_numeric(
                merged[column],
                errors="coerce"
            )

    feature_table["outcome"] = create_outcome(
        merged,
        survival_column
    )

    feature_table = feature_table.merge(
        gene_matrix,
        on="patient_id",
        how="left"
    )

    gene_feature_columns = [
        f"{gene}_mutated"
        for gene in selected_genes
    ]

    for column in gene_feature_columns:
        if column in feature_table.columns:
            feature_table[column] = (
                feature_table[column]
                .fillna(0)
                .astype(int)
            )

    feature_table = feature_table.drop_duplicates(
        subset=["patient_id"]
    )

    feature_table = feature_table[
        feature_table["outcome"].notna()
    ].copy()

    feature_table["outcome"] = (
        feature_table["outcome"]
        .astype(int)
    )

    return feature_table.reset_index(drop=True)


def create_feature_summary(
    feature_table: pd.DataFrame
) -> pd.DataFrame:
    """Tworzy raport jakości tabeli ML."""

    rows = []

    for column in feature_table.columns:
        rows.append(
            {
                "feature": column,
                "dtype": str(feature_table[column].dtype),
                "missing_values": int(
                    feature_table[column].isna().sum()
                ),
                "missing_percent": round(
                    feature_table[column].isna().mean() * 100,
                    2
                ),
                "unique_values": int(
                    feature_table[column].nunique(
                        dropna=True
                    )
                )
            }
        )

    return pd.DataFrame(rows)


# =====================================
# MAIN
# =====================================

def main() -> None:
    print("Loading merged clinical-genomic data...")

    merged = load_csv(MERGED_FILE)
    genomics = load_csv(GENOMICS_FILE)
    top_genes = load_csv(TOP_GENES_FILE)

    merged = standardize_patient_id(merged)
    genomics = standardize_patient_id(genomics)

    genomics["hugo_symbol"] = (
        genomics["hugo_symbol"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    print("Selecting genes for ML features...")

    selected_genes = choose_top_genes(
        top_genes=top_genes,
        genomics=genomics,
        top_n=TOP_N_GENES
    )

    print(
        "Selected genes: "
        + ", ".join(selected_genes)
    )

    print("Creating gene mutation matrix...")

    gene_matrix = create_gene_mutation_matrix(
        genomics,
        selected_genes
    )

    print("Building ML feature table...")

    feature_table = build_feature_table(
        merged,
        gene_matrix,
        selected_genes
    )

    feature_summary = create_feature_summary(
        feature_table
    )

    print("Saving outputs...")

    feature_table.to_csv(
        ML_OUTPUT_FILE,
        index=False
    )

    feature_summary.to_csv(
        ML_SUMMARY_FILE,
        index=False
    )

    print("\nFinished successfully.")
    print(
        f"Patients in ML table: "
        f"{len(feature_table):,}"
    )
    print(
        f"Features including ID and outcome: "
        f"{feature_table.shape[1]:,}"
    )

    print("\nOutcome distribution:")
    print(
        feature_table["outcome"]
        .value_counts()
        .sort_index()
    )

    print("\nSaved files:")
    print(ML_OUTPUT_FILE)
    print(ML_SUMMARY_FILE)


if __name__ == "__main__":
    main()