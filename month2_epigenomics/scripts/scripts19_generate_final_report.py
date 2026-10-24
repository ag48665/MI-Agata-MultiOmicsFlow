from pathlib import Path

import pandas as pd


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

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

FIGURES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "figures"
)

FINAL_REPORT_FILE = (
    REPORTS_DIR
    / "month2_final_report.md"
)


# =====================================
# INPUT FILES
# =====================================

INTEGRATION_SUMMARY_FILE = (
    TABLES_DIR
    / "epigenomic_integration_summary.csv"
)

MODEL_COMPARISON_FILE = (
    TABLES_DIR
    / "epigenomic_model_comparison.csv"
)

COX_METRICS_FILE = (
    TABLES_DIR
    / "cox_model_metrics.csv"
)

COX_COEFFICIENT_FILE = (
    TABLES_DIR
    / "cox_model_coefficients.csv"
)

FEATURE_IMPORTANCE_FILE = (
    TABLES_DIR
    / "combined_feature_importance.csv"
)

SELECTED_METHYLATION_FILE = (
    TABLES_DIR
    / "selected_methylation_features.csv"
)


# =====================================
# FUNCTIONS
# =====================================

def load_optional_csv(
    file_path: Path,
) -> pd.DataFrame | None:
    """Load CSV if it exists."""

    if not file_path.exists():
        print(
            f"Skipping missing file: "
            f"{file_path.name}"
        )
        return None

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def dataframe_to_markdown(
    df: pd.DataFrame,
    max_rows: int = 20,
) -> str:
    """Convert a dataframe to Markdown."""

    if df is None or df.empty:
        return "_No data available._"

    display_df = (
        df
        .head(max_rows)
        .copy()
    )

    try:
        return display_df.to_markdown(
            index=False
        )

    except ImportError:
        return display_df.to_string(
            index=False
        )


def create_figure_section() -> list[str]:
    """Create Markdown links to generated figures."""

    figure_names = [
        "beta_distribution.png",
        "methylation_pca.png",
        "clinical_genomic_confusion_matrix.png",
        "epigenomic_confusion_matrix.png",
        "clinical_genomic_feature_importance.png",
        "epigenomic_feature_importance.png",
        "cox_risk_groups_survival.png",
        "cox_top_hazard_ratios.png",
        "combined_feature_importance.png",
    ]

    lines = []

    for figure_name in figure_names:

        figure_path = (
            FIGURES_DIR
            / figure_name
        )

        if figure_path.exists():

            relative_path = (
                Path("..")
                / "figures"
                / figure_name
            )

            lines.extend(
                [
                    f"### {figure_name}",
                    "",
                    f"![{figure_name}]({relative_path.as_posix()})",
                    "",
                ]
            )

    return lines


def create_final_report() -> str:
    """Build final Month 2 report."""

    integration = load_optional_csv(
        INTEGRATION_SUMMARY_FILE
    )

    model_comparison = load_optional_csv(
        MODEL_COMPARISON_FILE
    )

    cox_metrics = load_optional_csv(
        COX_METRICS_FILE
    )

    cox_coefficients = load_optional_csv(
        COX_COEFFICIENT_FILE
    )

    feature_importance = load_optional_csv(
        FEATURE_IMPORTANCE_FILE
    )

    selected_methylation = load_optional_csv(
        SELECTED_METHYLATION_FILE
    )

    report_lines = [
        "# Month 2: Epigenomics Integration and AI-Ready Multi-Omics Pipeline Extension",
        "",
        "## Project overview",
        "",
        (
            "This module extends the Month 1 clinical-genomics workflow "
            "by incorporating DNA methylation data as an epigenomic layer."
        ),
        "",
        (
            "The workflow includes methylation data preparation, "
            "quality control, patient-level matching, multi-omics integration, "
            "feature selection, visualization, machine-learning modelling, "
            "survival analysis, and feature-importance interpretation."
        ),
        "",
        "## Workflow summary",
        "",
        "1. Methylation data preparation",
        "2. Methylation quality control",
        "3. TCGA sample-to-patient matching",
        "4. Clinical-genomic-epigenomic integration",
        "5. Epigenomic visualization",
        "6. High-variance methylation feature selection",
        "7. Random Forest baseline comparison",
        "8. Cox proportional hazards survival modelling",
        "9. Combined feature-importance analysis",
        "",
        "## Epigenomic integration",
        "",
        dataframe_to_markdown(
            integration
        ),
        "",
        "## Methylation feature selection",
        "",
        (
            "The methylation layer was reduced using variance-based "
            "feature selection to create an AI-ready feature matrix."
        ),
        "",
        "### Top selected methylation features",
        "",
        dataframe_to_markdown(
            selected_methylation,
            max_rows=20,
        ),
        "",
        "## Classification model comparison",
        "",
        (
            "A Random Forest classifier was used to compare the "
            "clinical-genomic baseline against the "
            "clinical-genomic-epigenomic feature set."
        ),
        "",
        dataframe_to_markdown(
            model_comparison
        ),
        "",
        "### Interpretation",
        "",
        (
            "The Month 2 model evaluates whether adding methylation "
            "features changes predictive performance relative to the "
            "Month 1 clinical-genomic baseline."
        ),
        "",
        (
            "Very high performance should be interpreted cautiously "
            "and assessed for overfitting or residual information leakage."
        ),
        "",
        "## Survival analysis",
        "",
        (
            "A penalized Cox proportional hazards model was trained "
            "using overall-survival duration and event information."
        ),
        "",
        "### Cox model metrics",
        "",
        dataframe_to_markdown(
            cox_metrics
        ),
        "",
        "### Top Cox model features",
        "",
        dataframe_to_markdown(
            cox_coefficients,
            max_rows=20,
        ),
        "",
        "## Feature importance",
        "",
        (
            "Feature-importance outputs from the classification and "
            "survival analyses were combined into a unified ranking."
        ),
        "",
        dataframe_to_markdown(
            feature_importance,
            max_rows=20,
        ),
        "",
        "## Generated figures",
        "",
    ]

    report_lines.extend(
        create_figure_section()
    )

    report_lines.extend(
        [
            "## Outputs generated in Month 2",
            "",
            "### Processed data",
            "",
            "- `data/processed/methylation/hm450_clean.csv`",
            "- `data/processed/methylation/methylation_patient_mapping.csv`",
            "- `data/processed/methylation/methylation_patient_matrix.csv`",
            "- `data/processed/clinical_genomic_epigenomic_merged.csv`",
            "- `data/processed/ml_epigenomic_feature_table.csv`",
            "",
            "### Tables",
            "",
            "- `outputs/tables/epigenomic_integration_summary.csv`",
            "- `outputs/tables/methylation_feature_variance.csv`",
            "- `outputs/tables/selected_methylation_features.csv`",
            "- `outputs/tables/epigenomic_model_comparison.csv`",
            "- `outputs/tables/cox_model_metrics.csv`",
            "- `outputs/tables/cox_model_coefficients.csv`",
            "- `outputs/tables/cox_patient_risk_scores.csv`",
            "- `outputs/tables/combined_feature_importance.csv`",
            "",
            "### Reports",
            "",
            "- `outputs/reports/methylation_data_overview.md`",
            "- `outputs/reports/methylation_qc_report.md`",
            "- `outputs/reports/methylation_patient_matching.md`",
            "- `outputs/reports/epigenomic_integration_report.md`",
            "- `outputs/reports/epigenomic_figures_report.md`",
            "- `outputs/reports/epigenomic_feature_selection.md`",
            "- `outputs/reports/epigenomic_baseline_model_report.md`",
            "- `outputs/reports/survival_model_report.md`",
            "- `outputs/reports/feature_importance_report.md`",
            "",
            "## Conclusion",
            "",
            (
                "Month 2 successfully adds DNA methylation to the existing "
                "clinical-genomic workflow and produces an integrated "
                "clinical-genomic-epigenomic dataset suitable for downstream "
                "machine-learning and survival analyses."
            ),
            "",
            (
                "The resulting pipeline is modular and provides the foundation "
                "for adding transcriptomic and additional omics layers in "
                "subsequent stages."
            ),
            "",
        ]
    )

    return "\n".join(
        report_lines
    )


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Generate final Month 2 report."""

    print(
        "Generating final Month 2 report..."
    )

    report = create_final_report()

    FINAL_REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print(
        "\nFinished successfully."
    )

    print(
        "\nFinal report saved to:\n"
        f"{FINAL_REPORT_FILE}"
    )


if __name__ == "__main__":
    main()