from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

TABLES_DIR = PROJECT_ROOT / "outputs" / "tables"
FIGURES_DIR = PROJECT_ROOT / "outputs" / "figures"
REPORTS_DIR = PROJECT_ROOT / "outputs" / "reports"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_TABLE = TABLES_DIR / "combined_feature_importance.csv"
OUTPUT_FIGURE = FIGURES_DIR / "combined_feature_importance.png"
REPORT_FILE = REPORTS_DIR / "feature_importance_report.md"

# =====================================================
# INPUT FILES
# =====================================================

FILES = [
    TABLES_DIR / "clinical_feature_importance.csv",
    TABLES_DIR / "genomic_feature_importance.csv",
    TABLES_DIR / "epigenomic_feature_importance.csv",
    TABLES_DIR / "cox_selected_features.csv",
]

# =====================================================
# FUNCTIONS
# =====================================================


def load_feature_tables():

    tables = []

    for file in FILES:

        if not file.exists():
            print(f"Skipping missing file: {file.name}")
            continue

        df = pd.read_csv(file)

        if "feature" not in df.columns:
            continue

        score_column = None

        for c in df.columns:
            if c.lower() in [
                "importance",
                "score",
                "coefficient",
                "coef",
                "abs_coef",
            ]:
                score_column = c
                break

        if score_column is None:

            numeric = df.select_dtypes(include="number").columns

            if len(numeric) == 0:
                continue

            score_column = numeric[0]

        tmp = pd.DataFrame(
            {
                "feature": df["feature"],
                "importance": df[score_column].abs(),
                "source": file.stem,
            }
        )

        tables.append(tmp)

    if len(tables) == 0:
        raise RuntimeError("No feature tables found.")

    return pd.concat(
        tables,
        ignore_index=True,
    )


def build_ranking(df):

    ranking = (
        df.groupby("feature")["importance"]
        .mean()
        .sort_values(ascending=False)
        .reset_index()
    )

    return ranking


def plot_importance(df):

    top = df.head(20)

    plt.figure(figsize=(9, 8))

    plt.barh(
        top["feature"][::-1],
        top["importance"][::-1],
    )

    plt.xlabel("Importance")

    plt.title("Top 20 most important features")

    plt.tight_layout()

    plt.savefig(
        OUTPUT_FIGURE,
        dpi=300,
    )

    plt.close()


def write_report(df):

    report = f"""
# Feature Importance Report

Total unique features: **{len(df):,}**

Top 10 features

{df.head(10).to_markdown(index=False)}
"""

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )


# =====================================================
# MAIN
# =====================================================


def main():

    print("Loading feature importance tables...")

    tables = load_feature_tables()

    print(f"Loaded rows: {len(tables):,}")

    ranking = build_ranking(tables)

    ranking.to_csv(
        OUTPUT_TABLE,
        index=False,
    )

    plot_importance(ranking)

    write_report(ranking)

    print("\nFinished successfully.")

    print(
        f"\nCombined table:\n{OUTPUT_TABLE}"
    )

    print(
        f"\nFigure:\n{OUTPUT_FIGURE}"
    )

    print(
        f"\nReport:\n{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()