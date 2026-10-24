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

METHYLATION_DIR = (
    PROCESSED_DIR
    / "methylation"
)

REPORTS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
)

REPORTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

METHYLATION_FILE = (
    METHYLATION_DIR
    / "hm450_clean.csv"
)

ML_FEATURE_FILE = (
    PROCESSED_DIR
    / "ml_feature_table.csv"
)

OUTPUT_FILE = (
    METHYLATION_DIR
    / "methylation_patient_mapping.csv"
)

REPORT_FILE = (
    REPORTS_DIR
    / "methylation_patient_matching.md"
)

# =====================================
# FUNCTIONS
# =====================================

def load_csv(file_path: Path) -> pd.DataFrame:

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def create_patient_id(sample_name: str) -> str:

    if pd.isna(sample_name):
        return ""

    return "-".join(
        str(sample_name).split("-")[:3]
    )

# =====================================
# MAIN
# =====================================

def main():

    print("Loading datasets...")

    methylation = load_csv(
        METHYLATION_FILE
    )

    ml = load_csv(
        ML_FEATURE_FILE
    )

    print(
        f"Methylation rows: {methylation.shape[0]:,}"
    )

    print(
        f"ML rows: {ml.shape[0]:,}"
    )

    sample_columns = [
        c
        for c in methylation.columns
        if c.startswith("TCGA-")
    ]

    patient_ids = [
        create_patient_id(c)
        for c in sample_columns
    ]

    mapping = pd.DataFrame({

        "sample_id": sample_columns,

        "patient_id": patient_ids

    })

    matched = mapping.merge(

        ml[["patient_id"]],

        on="patient_id",

        how="inner",

    )

    print(
        f"Matched patients: {matched.shape[0]}"
    )

    matched.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    report = f"""
# Methylation Patient Matching

Total methylation samples: {len(sample_columns)}

Matched patients: {matched.shape[0]}

Unmatched samples: {len(sample_columns)-matched.shape[0]}
"""

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print("\nFinished successfully.")

    print(
        f"Mapping saved to:\n{OUTPUT_FILE}"
    )

    print(
        f"Report saved to:\n{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()
