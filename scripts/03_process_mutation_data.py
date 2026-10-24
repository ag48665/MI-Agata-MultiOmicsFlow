from pathlib import Path

import pandas as pd


# =====================================
# ŚCIEŻKI
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

MUTATION_FILE = RAW_DATA_DIR / "data_mutations.txt"


# =====================================
# FUNKCJE
# =====================================

def load_mutation_data(file_path: Path) -> pd.DataFrame:
    """Wczytuje dane mutacji z pliku cBioPortal."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku mutacji:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        low_memory=False
    )


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standaryzuje nazwy kolumn."""

    df = df.copy()

    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
        .str.replace("-", "_", regex=False)
    )

    return df


def clean_gene_names(df: pd.DataFrame) -> pd.DataFrame:
    """Czyści i standaryzuje symbole genów."""

    df = df.copy()

    if "hugo_symbol" not in df.columns:
        raise KeyError(
            "Nie znaleziono kolumny 'hugo_symbol'."
        )

    df["hugo_symbol"] = (
        df["hugo_symbol"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    df = df[
        df["hugo_symbol"].notna()
        & (df["hugo_symbol"] != "")
        & (df["hugo_symbol"] != "NAN")
    ].copy()

    return df


def standardize_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Tworzy ujednolicone identyfikatory próbki i pacjenta
    na podstawie kodu TCGA.
    """

    df = df.copy()

    if "tumor_sample_barcode" not in df.columns:
        raise KeyError(
            "Nie znaleziono kolumny 'tumor_sample_barcode'."
        )

    df["sample_id"] = (
        df["tumor_sample_barcode"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    df["patient_id"] = (
        df["sample_id"]
        .str.split("-")
        .str[:3]
        .str.join("-")
    )

    return df


def clean_variant_classification(df: pd.DataFrame) -> pd.DataFrame:
    """Czyści nazwy klasyfikacji wariantów."""

    df = df.copy()

    if "variant_classification" not in df.columns:
        raise KeyError(
            "Nie znaleziono kolumny 'variant_classification'."
        )

    df["variant_classification"] = (
        df["variant_classification"]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
    )

    return df


def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Usuwa całkowicie zduplikowane rekordy."""

    before = len(df)
    df = df.drop_duplicates().copy()
    removed = before - len(df)

    print(f"Removed duplicate records: {removed:,}")

    return df


def mutation_counts_per_patient(df: pd.DataFrame) -> pd.DataFrame:
    """Zlicza liczbę rekordów mutacji na pacjenta."""

    return (
        df.groupby("patient_id", dropna=False)
        .size()
        .reset_index(name="mutation_count")
        .sort_values(
            "mutation_count",
            ascending=False
        )
        .reset_index(drop=True)
    )


def top_mutated_genes(
    df: pd.DataFrame,
    top_n: int = 20
) -> pd.DataFrame:
    """Zwraca najczęściej mutowane geny."""

    summary = (
        df.groupby("hugo_symbol")
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

    return summary


def variant_classification_summary(
    df: pd.DataFrame
) -> pd.DataFrame:
    """Tworzy podsumowanie typów mutacji."""

    return (
        df.groupby(
            "variant_classification",
            dropna=False
        )
        .size()
        .reset_index(name="count")
        .sort_values(
            "count",
            ascending=False
        )
        .reset_index(drop=True)
    )


def select_cleaned_columns(
    df: pd.DataFrame
) -> pd.DataFrame:
    """
    Wybiera najważniejsze kolumny potrzebne
    w kolejnych etapach pipeline'u.
    """

    preferred_columns = [
        "patient_id",
        "sample_id",
        "hugo_symbol",
        "entrez_gene_id",
        "ncbi_build",
        "chromosome",
        "start_position",
        "end_position",
        "strand",
        "variant_classification",
        "variant_type",
        "reference_allele",
        "tumor_seq_allele1",
        "tumor_seq_allele2",
        "tumor_sample_barcode",
        "matched_norm_sample_barcode",
        "hgvsc",
        "hgvsp",
        "hgvsp_short",
        "protein_position",
        "codons",
        "t_ref_count",
        "t_alt_count",
        "n_ref_count",
        "n_alt_count",
        "annotation_status"
    ]

    available_columns = [
        column
        for column in preferred_columns
        if column in df.columns
    ]

    return df[available_columns].copy()


# =====================================
# MAIN
# =====================================

def main() -> None:
    print("Loading mutation data...")

    mutations = load_mutation_data(MUTATION_FILE)

    print(
        f"Loaded {len(mutations):,} records "
        f"and {mutations.shape[1]} columns."
    )

    print("Standardizing columns...")
    mutations = standardize_columns(mutations)

    print("Cleaning gene names...")
    mutations = clean_gene_names(mutations)

    print("Standardizing TCGA identifiers...")
    mutations = standardize_ids(mutations)

    print("Cleaning variant classifications...")
    mutations = clean_variant_classification(mutations)

    print("Removing duplicates...")
    mutations = remove_duplicates(mutations)

    print("Generating summaries...")

    mutation_counts = mutation_counts_per_patient(
        mutations
    )

    top_genes = top_mutated_genes(
        mutations,
        top_n=20
    )

    variant_summary = variant_classification_summary(
        mutations
    )

    genomics_cleaned = select_cleaned_columns(
        mutations
    )

    print("Saving outputs...")

    genomics_cleaned.to_csv(
        PROCESSED_DIR / "genomics_cleaned.csv",
        index=False
    )

    mutation_counts.to_csv(
        TABLES_DIR / "mutation_counts_per_patient.csv",
        index=False
    )

    top_genes.to_csv(
        TABLES_DIR / "top_mutated_genes.csv",
        index=False
    )

    variant_summary.to_csv(
        TABLES_DIR / "variant_classification_summary.csv",
        index=False
    )

    print("\nFinished successfully.")
    print(f"Cleaned mutation records: {len(genomics_cleaned):,}")
    print(
        f"Unique patients: "
        f"{genomics_cleaned['patient_id'].nunique():,}"
    )
    print(
        f"Unique samples: "
        f"{genomics_cleaned['sample_id'].nunique():,}"
    )
    print(
        f"Unique genes: "
        f"{genomics_cleaned['hugo_symbol'].nunique():,}"
    )

    print("\nSaved files:")
    print(PROCESSED_DIR / "genomics_cleaned.csv")
    print(TABLES_DIR / "mutation_counts_per_patient.csv")
    print(TABLES_DIR / "top_mutated_genes.csv")
    print(TABLES_DIR / "variant_classification_summary.csv")


if __name__ == "__main__":
    main()