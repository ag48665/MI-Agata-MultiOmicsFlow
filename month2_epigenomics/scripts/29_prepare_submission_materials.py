from pathlib import Path
import hashlib
import shutil

import pandas as pd


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUTS_DIR = (
    PROJECT_ROOT
    / "outputs"
)

TABLES_DIR = (
    OUTPUTS_DIR
    / "tables"
)

FIGURES_DIR = (
    OUTPUTS_DIR
    / "figures"
)

REPORTS_DIR = (
    OUTPUTS_DIR
    / "reports"
)

SUBMISSION_DIR = (
    OUTPUTS_DIR
    / "submission"
)

MAIN_FIGURES_DIR = (
    SUBMISSION_DIR
    / "main_figures"
)

MAIN_TABLES_DIR = (
    SUBMISSION_DIR
    / "main_tables"
)

SUPPLEMENTARY_DIR = (
    SUBMISSION_DIR
    / "supplementary"
)

SUPPLEMENTARY_FIGURES_DIR = (
    SUPPLEMENTARY_DIR
    / "figures"
)

SUPPLEMENTARY_TABLES_DIR = (
    SUPPLEMENTARY_DIR
    / "tables"
)


MANIFEST_FILE = (
    SUBMISSION_DIR
    / "submission_manifest.csv"
)

FIGURE_CAPTIONS_FILE = (
    SUBMISSION_DIR
    / "figure_captions.md"
)

TABLE_CAPTIONS_FILE = (
    SUBMISSION_DIR
    / "table_captions.md"
)

MANUSCRIPT_OUTLINE_FILE = (
    SUBMISSION_DIR
    / "manuscript_outline.md"
)

README_FILE = (
    SUBMISSION_DIR
    / "README.md"
)


# =====================================================
# MAIN FIGURES
# =====================================================

MAIN_FIGURES = {
    "Figure_1_model_validation.png":
        "model_validation_roc_auc.png",

    "Figure_2_precision_recall.png":
        "threshold_validation_precision_recall.png",

    "Figure_3_calibration.png":
        "model_calibration_curve.png",

    "Figure_4_confusion_matrix.png":
        "confusion_matrix_optimized_threshold.png",

    "Figure_5_survival.png":
        "cox_risk_groups_survival.png",

    "Figure_6_feature_importance.png":
        "combined_feature_importance.png",
}


# =====================================================
# MAIN TABLES
# =====================================================

MAIN_TABLES = {
    "Table_1_model_validation.csv":
        "model_validation_summary.csv",

    "Table_2_threshold_optimization.csv":
        "threshold_optimization_summary.csv",

    "Table_3_calibration.csv":
        "calibration_validation_summary.csv",

    "Table_4_final_model_metrics.csv":
        "final_epigenomic_model_metrics.csv",

    "Table_5_key_results.csv":
        "month2_key_results.csv",
}


# =====================================================
# SUPPLEMENTARY FIGURES
# =====================================================

SUPPLEMENTARY_FIGURES = [
    "beta_distribution.png",
    "methylation_pca.png",
    "threshold_validation_roc_curve.png",
    "confusion_matrix_threshold_050.png",
    "calibrated_probability_distribution.png",
    "cox_top_hazard_ratios.png",
]


# =====================================================
# SUPPLEMENTARY TABLES
# =====================================================

SUPPLEMENTARY_TABLES = [
    "model_validation_folds.csv",
    "model_leakage_audit.csv",
    "cv_selected_methylation_features.csv",
    "optimized_thresholds.csv",
    "threshold_validation_folds.csv",
    "threshold_validation_predictions.csv",
    "calibration_validation_folds.csv",
    "calibration_predictions.csv",
    "final_calibrated_oof_predictions.csv",
    "final_calibrated_threshold_search.csv",
    "final_epigenomic_model_features.csv",
    "month2_final_qc.csv",
]


# =====================================================
# HELPERS
# =====================================================

def reset_submission_directory() -> None:
    """Rebuild submission directory."""

    if SUBMISSION_DIR.exists():
        shutil.rmtree(
            SUBMISSION_DIR
        )

    for directory in [
        SUBMISSION_DIR,
        MAIN_FIGURES_DIR,
        MAIN_TABLES_DIR,
        SUPPLEMENTARY_DIR,
        SUPPLEMENTARY_FIGURES_DIR,
        SUPPLEMENTARY_TABLES_DIR,
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def copy_file(
    source: Path,
    destination: Path,
) -> bool:
    """Copy file if it exists."""

    if not source.exists():
        print(
            f"Missing, skipped: {source.name}"
        )
        return False

    shutil.copy2(
        source,
        destination,
    )

    print(
        f"Copied: {destination.name}"
    )

    return True


def sha256_file(
    file_path: Path,
) -> str:
    """Calculate SHA256 checksum."""

    digest = hashlib.sha256()

    with file_path.open(
        "rb"
    ) as handle:

        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(
                chunk
            )

    return digest.hexdigest()


# =====================================================
# COPY MAIN MATERIALS
# =====================================================

def copy_main_figures() -> int:
    """Copy publication main figures."""

    count = 0

    for new_name, source_name in (
        MAIN_FIGURES.items()
    ):

        source = (
            FIGURES_DIR
            / source_name
        )

        destination = (
            MAIN_FIGURES_DIR
            / new_name
        )

        if copy_file(
            source,
            destination,
        ):
            count += 1

    return count


def copy_main_tables() -> int:
    """Copy publication main tables."""

    count = 0

    for new_name, source_name in (
        MAIN_TABLES.items()
    ):

        source = (
            TABLES_DIR
            / source_name
        )

        destination = (
            MAIN_TABLES_DIR
            / new_name
        )

        if copy_file(
            source,
            destination,
        ):
            count += 1

    return count


# =====================================================
# COPY SUPPLEMENTARY MATERIALS
# =====================================================

def copy_supplementary_figures() -> int:
    """Copy supplementary figures."""

    count = 0

    for filename in (
        SUPPLEMENTARY_FIGURES
    ):

        source = (
            FIGURES_DIR
            / filename
        )

        destination = (
            SUPPLEMENTARY_FIGURES_DIR
            / filename
        )

        if copy_file(
            source,
            destination,
        ):
            count += 1

    return count


def copy_supplementary_tables() -> int:
    """Copy supplementary tables."""

    count = 0

    for filename in (
        SUPPLEMENTARY_TABLES
    ):

        source = (
            TABLES_DIR
            / filename
        )

        destination = (
            SUPPLEMENTARY_TABLES_DIR
            / filename
        )

        if copy_file(
            source,
            destination,
        ):
            count += 1

    return count


# =====================================================
# CAPTIONS
# =====================================================

def create_figure_captions() -> str:
    """Create manuscript-style figure captions."""

    lines = [
        "# Figure Captions",
        "",
        "## Figure 1",
        "",
        (
            "**Leakage-aware model validation.** "
            "Comparison of clinical-genomic and "
            "clinical-genomic-epigenomic model performance "
            "across cross-validation folds. Inclusion of DNA "
            "methylation features increased discriminatory "
            "performance relative to the baseline model."
        ),
        "",
        "## Figure 2",
        "",
        (
            "**Precision-recall performance of the "
            "epigenomic classifier.** "
            "Precision-recall curve derived from leakage-aware "
            "out-of-fold probability estimates. This figure "
            "illustrates model performance under class imbalance."
        ),
        "",
        "## Figure 3",
        "",
        (
            "**Probability calibration.** "
            "Calibration curves comparing uncalibrated Random "
            "Forest probabilities with sigmoid and isotonic "
            "calibration. The diagonal represents perfect "
            "agreement between predicted probabilities and "
            "observed outcome frequencies."
        ),
        "",
        "## Figure 4",
        "",
        (
            "**Confusion matrix at the optimized operating "
            "threshold.** "
            "Classification results after threshold optimisation "
            "designed to improve balanced accuracy and minority-"
            "class sensitivity."
        ),
        "",
        "## Figure 5",
        "",
        (
            "**Kaplan-Meier survival analysis by model-derived "
            "risk group.** "
            "Patients are stratified according to model-derived "
            "risk estimates, illustrating differences in survival "
            "between predicted risk groups."
        ),
        "",
        "## Figure 6",
        "",
        (
            "**Combined feature importance.** "
            "Summary of the most influential features identified "
            "across the multi-omics modelling workflow."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


def create_table_captions() -> str:
    """Create manuscript-style table captions."""

    lines = [
        "# Table Captions",
        "",
        "## Table 1",
        "",
        (
            "**Leakage-aware cross-validation performance.** "
            "Comparison of clinical-genomic and "
            "clinical-genomic-epigenomic models across five "
            "stratified cross-validation folds."
        ),
        "",
        "## Table 2",
        "",
        (
            "**Classification-threshold optimisation.** "
            "Comparison of performance at the default threshold "
            "of 0.50 and the optimized training-derived threshold."
        ),
        "",
        "## Table 3",
        "",
        (
            "**Probability calibration performance.** "
            "Comparison of uncalibrated, sigmoid-calibrated and "
            "isotonic-calibrated probability estimates using "
            "ROC AUC, PR AUC, Brier score and log loss."
        ),
        "",
        "## Table 4",
        "",
        (
            "**Final calibrated model performance.** "
            "Leakage-aware out-of-fold metrics for the final "
            "Random Forest model with isotonic probability "
            "calibration."
        ),
        "",
        "## Table 5",
        "",
        (
            "**Summary of principal Month 2 results.** "
            "Key model-validation, threshold-optimisation and "
            "calibration results used in the manuscript."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


# =====================================================
# MANUSCRIPT OUTLINE
# =====================================================

def create_manuscript_outline() -> str:
    """Create recommended manuscript structure."""

    lines = [
        "# Manuscript Outline",
        "",
        "## 1. Introduction",
        "",
        (
            "Introduce breast cancer heterogeneity, the value "
            "of multi-omics integration and the potential "
            "predictive contribution of DNA methylation."
        ),
        "",
        "## 2. Materials and Methods",
        "",
        "### 2.1 Cohort and data sources",
        "### 2.2 Clinical data preprocessing",
        "### 2.3 Genomic data preprocessing",
        "### 2.4 DNA methylation preprocessing and QC",
        "### 2.5 Patient-level multi-omics integration",
        "### 2.6 Feature selection",
        "### 2.7 Machine-learning model development",
        "### 2.8 Leakage-aware cross-validation",
        "### 2.9 Threshold optimisation",
        "### 2.10 Probability calibration",
        "### 2.11 Survival analysis",
        "### 2.12 Statistical and computational implementation",
        "",
        "## 3. Results",
        "",
        "### 3.1 Cohort and integrated dataset",
        "### 3.2 Contribution of methylation to prediction",
        "### 3.3 Threshold optimisation",
        "### 3.4 Probability calibration",
        "### 3.5 Final model performance",
        "### 3.6 Survival analysis",
        "### 3.7 Feature importance",
        "",
        (
            "Primary source text: "
            "`outputs/reports/month2_manuscript_results.md`"
        ),
        "",
        "## 4. Discussion",
        "",
        (
            "Primary source text: "
            "`outputs/reports/month2_manuscript_discussion.md`"
        ),
        "",
        "Discuss:",
        "",
        "- added value of methylation;",
        "- reduction of performance after leakage-aware validation;",
        "- importance of class imbalance;",
        "- threshold trade-offs;",
        "- probability calibration;",
        "- limitations of internal validation;",
        "- need for external validation.",
        "",
        "## 5. Conclusions",
        "",
        (
            "Primary source text: "
            "`outputs/reports/month2_manuscript_conclusions.md`"
        ),
        "",
        "## Main figures",
        "",
        "1. Model-validation comparison",
        "2. Precision-recall curve",
        "3. Calibration curve",
        "4. Optimized confusion matrix",
        "5. Kaplan-Meier survival figure",
        "6. Feature-importance figure",
        "",
        "## Main tables",
        "",
        "1. Model-validation metrics",
        "2. Threshold comparison",
        "3. Calibration comparison",
        "4. Final-model metrics",
        "5. Key results",
        "",
        "## Supplementary material",
        "",
        (
            "Include detailed fold-level results, leakage audit, "
            "selected methylation features, OOF predictions, "
            "threshold-search results and QC tables."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


# =====================================================
# MANIFEST
# =====================================================

def create_manifest() -> pd.DataFrame:
    """Create submission-file manifest."""

    rows = []

    for file_path in sorted(
        SUBMISSION_DIR.rglob("*")
    ):

        if not file_path.is_file():
            continue

        if file_path.name == (
            MANIFEST_FILE.name
        ):
            continue

        relative = (
            file_path.relative_to(
                SUBMISSION_DIR
            )
        )

        rows.append(
            {
                "file":
                    relative.as_posix(),
                "size_bytes":
                    file_path.stat().st_size,
                "sha256":
                    sha256_file(
                        file_path
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =====================================================
# README
# =====================================================

def create_readme(
    main_figures_count: int,
    main_tables_count: int,
    supplementary_figures_count: int,
    supplementary_tables_count: int,
) -> str:
    """Create submission package README."""

    lines = [
        "# Submission Materials",
        "",
        (
            "Prepared publication-oriented materials from "
            "MI-Agata-MultiOmicsFlow Month 2."
        ),
        "",
        "## Main manuscript materials",
        "",
        (
            f"- Main figures: "
            f"{main_figures_count}"
        ),
        (
            f"- Main tables: "
            f"{main_tables_count}"
        ),
        "",
        "## Supplementary materials",
        "",
        (
            f"- Supplementary figures: "
            f"{supplementary_figures_count}"
        ),
        (
            f"- Supplementary tables: "
            f"{supplementary_tables_count}"
        ),
        "",
        "## Supporting files",
        "",
        "- `figure_captions.md`",
        "- `table_captions.md`",
        "- `manuscript_outline.md`",
        "- `submission_manifest.csv`",
        "",
        "## Important note",
        "",
        (
            "These materials summarize internally validated "
            "development-stage results. They do not represent "
            "external clinical validation."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Prepare publication/submission materials."""

    print(
        "Preparing submission materials..."
    )

    reset_submission_directory()

    print(
        "\nCopying main figures..."
    )

    main_figures_count = (
        copy_main_figures()
    )

    print(
        "\nCopying main tables..."
    )

    main_tables_count = (
        copy_main_tables()
    )

    print(
        "\nCopying supplementary figures..."
    )

    supplementary_figures_count = (
        copy_supplementary_figures()
    )

    print(
        "\nCopying supplementary tables..."
    )

    supplementary_tables_count = (
        copy_supplementary_tables()
    )

    # =================================================
    # WRITE CAPTIONS
    # =================================================

    FIGURE_CAPTIONS_FILE.write_text(
        create_figure_captions(),
        encoding="utf-8",
    )

    TABLE_CAPTIONS_FILE.write_text(
        create_table_captions(),
        encoding="utf-8",
    )

    MANUSCRIPT_OUTLINE_FILE.write_text(
        create_manuscript_outline(),
        encoding="utf-8",
    )

    # =================================================
    # README
    # =================================================

    README_FILE.write_text(
        create_readme(
            main_figures_count=(
                main_figures_count
            ),
            main_tables_count=(
                main_tables_count
            ),
            supplementary_figures_count=(
                supplementary_figures_count
            ),
            supplementary_tables_count=(
                supplementary_tables_count
            ),
        ),
        encoding="utf-8",
    )

    # =================================================
    # MANIFEST
    # =================================================

    manifest = create_manifest()

    manifest.to_csv(
        MANIFEST_FILE,
        index=False,
    )

    # =================================================
    # DONE
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "SUBMISSION MATERIALS COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        f"\nMain figures: "
        f"{main_figures_count}"
    )

    print(
        f"Main tables: "
        f"{main_tables_count}"
    )

    print(
        f"Supplementary figures: "
        f"{supplementary_figures_count}"
    )

    print(
        f"Supplementary tables: "
        f"{supplementary_tables_count}"
    )

    print(
        "\nSubmission folder saved to:\n"
        f"{SUBMISSION_DIR}"
    )

    print(
        "\nFigure captions:\n"
        f"{FIGURE_CAPTIONS_FILE}"
    )

    print(
        "\nTable captions:\n"
        f"{TABLE_CAPTIONS_FILE}"
    )

    print(
        "\nManuscript outline:\n"
        f"{MANUSCRIPT_OUTLINE_FILE}"
    )

    print(
        "\nManifest:\n"
        f"{MANIFEST_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()