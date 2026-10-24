from pathlib import Path
import hashlib
import json
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

MODELS_DIR = (
    OUTPUTS_DIR
    / "models"
)

PACKAGE_DIR = (
    OUTPUTS_DIR
    / "month2_package"
)

PACKAGE_TABLES_DIR = (
    PACKAGE_DIR
    / "tables"
)

PACKAGE_FIGURES_DIR = (
    PACKAGE_DIR
    / "figures"
)

PACKAGE_REPORTS_DIR = (
    PACKAGE_DIR
    / "reports"
)

PACKAGE_MODELS_DIR = (
    PACKAGE_DIR
    / "models"
)

FINAL_REPORT_FILE = (
    REPORTS_DIR
    / "month2_final_report.md"
)

PACKAGE_README_FILE = (
    PACKAGE_DIR
    / "README.md"
)

PACKAGE_MANIFEST_FILE = (
    PACKAGE_DIR
    / "manifest.csv"
)

PACKAGE_METADATA_FILE = (
    PACKAGE_DIR
    / "package_metadata.json"
)


# =====================================================
# IMPORTANT RESULT FILES
# =====================================================

VALIDATION_SUMMARY_FILE = (
    TABLES_DIR
    / "model_validation_summary.csv"
)

THRESHOLD_SUMMARY_FILE = (
    TABLES_DIR
    / "threshold_optimization_summary.csv"
)

CALIBRATION_SUMMARY_FILE = (
    TABLES_DIR
    / "calibration_validation_summary.csv"
)

FINAL_METRICS_FILE = (
    TABLES_DIR
    / "final_epigenomic_model_metrics.csv"
)

FINAL_THRESHOLD_FILE = (
    TABLES_DIR
    / "final_calibrated_threshold.csv"
)

FINAL_FEATURES_FILE = (
    TABLES_DIR
    / "final_epigenomic_model_features.csv"
)

FINAL_MODEL_CONFIG_FILE = (
    MODELS_DIR
    / "final_epigenomic_model_config.json"
)


# =====================================================
# HELPERS
# =====================================================

def load_optional_csv(
    file_path: Path,
) -> pd.DataFrame | None:
    """Load CSV if available."""

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
    df: pd.DataFrame | None,
    max_rows: int = 20,
) -> str:
    """Convert dataframe to Markdown safely."""

    if df is None or df.empty:
        return "_Data not available._"

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
        return (
            "```text\n"
            + display_df.to_string(
                index=False
            )
            + "\n```"
        )


def metric_from_table(
    df: pd.DataFrame | None,
    metric_name: str,
) -> float | None:
    """Read metric from metric/value table."""

    if df is None:
        return None

    if (
        "metric" not in df.columns
        or "value" not in df.columns
    ):
        return None

    rows = df[
        df["metric"]
        == metric_name
    ]

    if rows.empty:
        return None

    try:
        return float(
            rows.iloc[0]["value"]
        )

    except (TypeError, ValueError):
        return None


def format_metric(
    value: float | None,
    digits: int = 3,
) -> str:
    """Format optional metric."""

    if value is None:
        return "N/A"

    return f"{value:.{digits}f}"


# =====================================================
# FINAL REPORT
# =====================================================

def create_updated_final_report() -> str:
    """Create final Month 2 report including Steps 21-24."""

    validation = load_optional_csv(
        VALIDATION_SUMMARY_FILE
    )

    threshold = load_optional_csv(
        THRESHOLD_SUMMARY_FILE
    )

    calibration = load_optional_csv(
        CALIBRATION_SUMMARY_FILE
    )

    final_metrics = load_optional_csv(
        FINAL_METRICS_FILE
    )

    final_threshold = load_optional_csv(
        FINAL_THRESHOLD_FILE
    )

    final_features = load_optional_csv(
        FINAL_FEATURES_FILE
    )

    roc_auc = metric_from_table(
        final_metrics,
        "roc_auc",
    )

    pr_auc = metric_from_table(
        final_metrics,
        "pr_auc",
    )

    brier = metric_from_table(
        final_metrics,
        "brier_score",
    )

    balanced_accuracy = metric_from_table(
        final_metrics,
        "balanced_accuracy",
    )

    precision = metric_from_table(
        final_metrics,
        "precision",
    )

    recall = metric_from_table(
        final_metrics,
        "recall",
    )

    f1 = metric_from_table(
        final_metrics,
        "f1_score",
    )

    threshold_value = metric_from_table(
        final_metrics,
        "threshold",
    )

    report_lines = [
        "# Month 2 Final Report",
        "",
        "## Epigenomics Integration and "
        "AI-Ready Multi-Omics Pipeline Extension",
        "",
        "## Overview",
        "",
        (
            "Month 2 extended the clinical-genomic "
            "workflow with DNA methylation data."
        ),
        "",
        (
            "The final workflow includes methylation "
            "quality control, patient matching, "
            "multi-omics integration, feature selection, "
            "machine-learning classification, survival "
            "modelling, leakage-aware validation, "
            "decision-threshold optimisation, probability "
            "calibration and final model packaging."
        ),
        "",
        "## Completed workflow",
        "",
        "1. Methylation data preparation",
        "2. Methylation quality control",
        "3. Sample-to-patient matching",
        "4. Clinical-genomic-epigenomic integration",
        "5. Epigenomic visualisation",
        "6. Methylation feature selection",
        "7. AI-ready feature-table generation",
        "8. Classification modelling",
        "9. Cox survival modelling",
        "10. Feature-importance analysis",
        "11. Leakage-aware model validation",
        "12. Classification-threshold optimisation",
        "13. Probability calibration",
        "14. Final calibrated model training",
        "",
        "## Leakage-aware model validation",
        "",
        (
            "Step 21 re-evaluated classification performance "
            "using stratified cross-validation and leakage "
            "controls."
        ),
        "",
        dataframe_to_markdown(
            validation
        ),
        "",
        "### Interpretation",
        "",
        (
            "The cross-validation analysis provides a more "
            "reliable estimate of generalisation performance "
            "than the original single train/test result."
        ),
        "",
        (
            "Patient identifiers, sample identifiers, survival "
            "variables and obvious outcome proxies were excluded "
            "from the predictor set."
        ),
        "",
        "## Decision-threshold optimisation",
        "",
        (
            "Because the outcome classes are imbalanced, "
            "classification performance was evaluated at both "
            "the conventional threshold of 0.50 and a threshold "
            "selected from training-fold predictions."
        ),
        "",
        dataframe_to_markdown(
            threshold
        ),
        "",
        (
            "Threshold optimisation substantially changed the "
            "trade-off between precision and recall. "
            "ROC AUC and PR AUC remained threshold-independent."
        ),
        "",
        "## Probability calibration",
        "",
        (
            "Uncalibrated Random Forest probabilities were "
            "compared with sigmoid and isotonic calibration."
        ),
        "",
        dataframe_to_markdown(
            calibration
        ),
        "",
        (
            "Isotonic calibration achieved the lowest mean "
            "Brier score during cross-validation and was "
            "therefore selected for the final development model."
        ),
        "",
        "## Final calibrated model",
        "",
        (
            "The final development-stage classifier is a "
            "Random Forest model with isotonic probability "
            "calibration."
        ),
        "",
        (
            "The final operating threshold was estimated "
            "using leakage-aware out-of-fold calibrated "
            "probabilities."
        ),
        "",
        "### Final internal performance",
        "",
        f"- ROC AUC: **{format_metric(roc_auc)}**",
        f"- PR AUC: **{format_metric(pr_auc)}**",
        (
            "- Brier score: "
            f"**{format_metric(brier, 4)}**"
        ),
        (
            "- Balanced accuracy: "
            f"**{format_metric(balanced_accuracy)}**"
        ),
        (
            "- Precision: "
            f"**{format_metric(precision)}**"
        ),
        (
            "- Recall: "
            f"**{format_metric(recall)}**"
        ),
        (
            "- F1-score: "
            f"**{format_metric(f1)}**"
        ),
        (
            "- Final calibrated threshold: "
            f"**{format_metric(threshold_value)}**"
        ),
        "",
        "### Final feature set",
        "",
        dataframe_to_markdown(
            final_features,
            max_rows=25,
        ),
        "",
        "## Final model files",
        "",
        "- `outputs/models/final_epigenomic_calibrated_model.pkl`",
        "- `outputs/models/final_epigenomic_model_config.json`",
        "",
        "## Final validation tables",
        "",
        "- `outputs/tables/model_validation_summary.csv`",
        "- `outputs/tables/model_validation_folds.csv`",
        "- `outputs/tables/model_leakage_audit.csv`",
        "- `outputs/tables/threshold_optimization_summary.csv`",
        "- `outputs/tables/optimized_thresholds.csv`",
        "- `outputs/tables/calibration_validation_summary.csv`",
        "- `outputs/tables/calibration_validation_folds.csv`",
        "- `outputs/tables/final_epigenomic_model_metrics.csv`",
        "- `outputs/tables/final_calibrated_threshold.csv`",
        "- `outputs/tables/final_calibrated_oof_predictions.csv`",
        "",
        "## Final figures",
        "",
        "![Validation ROC AUC](../figures/model_validation_roc_auc.png)",
        "",
        "![Threshold ROC curve](../figures/threshold_validation_roc_curve.png)",
        "",
        "![Precision-recall curve](../figures/threshold_validation_precision_recall.png)",
        "",
        "![Calibration curve](../figures/model_calibration_curve.png)",
        "",
        "![Calibrated probability distribution](../figures/calibrated_probability_distribution.png)",
        "",
        "![Optimized confusion matrix](../figures/confusion_matrix_optimized_threshold.png)",
        "",
        "## Methodological interpretation",
        "",
        (
            "The final internal results support the conclusion "
            "that DNA methylation contributes additional "
            "predictive information beyond the original "
            "clinical-genomic feature set."
        ),
        "",
        (
            "However, all reported performance values remain "
            "internal development estimates."
        ),
        "",
        (
            "The final calibrated probabilities and operating "
            "threshold must not be interpreted as externally "
            "validated clinical risk estimates."
        ),
        "",
        (
            "Independent external validation on a separate "
            "cohort is required before clinical use."
        ),
        "",
        "## Conclusion",
        "",
        (
            "Month 2 produced a complete clinical-genomic-"
            "epigenomic machine-learning workflow with explicit "
            "controls for data leakage, class imbalance, "
            "threshold selection and probability calibration."
        ),
        "",
        (
            "The resulting calibrated model and associated "
            "validation artefacts provide a reproducible "
            "development-stage foundation for future external "
            "validation and additional multi-omics integration."
        ),
        "",
    ]

    return "\n".join(
        report_lines
    )


# =====================================================
# PACKAGE COPYING
# =====================================================

def prepare_package() -> None:
    """Rebuild package directory."""

    if PACKAGE_DIR.exists():
        shutil.rmtree(
            PACKAGE_DIR
        )

    for directory in [
        PACKAGE_DIR,
        PACKAGE_TABLES_DIR,
        PACKAGE_FIGURES_DIR,
        PACKAGE_REPORTS_DIR,
        PACKAGE_MODELS_DIR,
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def copy_all_files(
    source_dir: Path,
    destination_dir: Path,
    extensions: tuple[str, ...],
) -> int:
    """Copy all matching files."""

    count = 0

    if not source_dir.exists():
        return count

    for file_path in sorted(
        source_dir.iterdir()
    ):

        if not file_path.is_file():
            continue

        if (
            file_path.suffix.lower()
            not in extensions
        ):
            continue

        shutil.copy2(
            file_path,
            destination_dir
            / file_path.name,
        )

        count += 1

    return count


# =====================================================
# MANIFEST
# =====================================================

def calculate_sha256(
    file_path: Path,
) -> str:
    """Calculate SHA256."""

    digest = hashlib.sha256()

    with file_path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def create_manifest() -> pd.DataFrame:
    """Create package manifest."""

    rows = []

    for file_path in sorted(
        PACKAGE_DIR.rglob("*")
    ):

        if not file_path.is_file():
            continue

        if file_path.name in {
            "manifest.csv",
            "package_metadata.json",
        }:
            continue

        relative = file_path.relative_to(
            PACKAGE_DIR
        )

        rows.append(
            {
                "file":
                    relative.as_posix(),
                "size_bytes":
                    file_path.stat().st_size,
                "sha256":
                    calculate_sha256(
                        file_path
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =====================================================
# PACKAGE README
# =====================================================

def create_package_readme(
    reports_count: int,
    tables_count: int,
    figures_count: int,
    models_count: int,
) -> str:
    """Create final package README."""

    lines = [
        "# MI-Agata-MultiOmicsFlow",
        "",
        "## Month 2 Final Results Package",
        "",
        (
            "This package contains the complete Month 2 "
            "epigenomics extension of the multi-omics "
            "machine-learning workflow."
        ),
        "",
        "## Package contents",
        "",
        (
            f"- Reports: {reports_count}"
        ),
        (
            f"- Tables: {tables_count}"
        ),
        (
            f"- Figures: {figures_count}"
        ),
        (
            f"- Models/configuration files: "
            f"{models_count}"
        ),
        "",
        "## Main report",
        "",
        "`reports/month2_final_report.md`",
        "",
        "## Final model",
        "",
        (
            "`models/"
            "final_epigenomic_calibrated_model.pkl`"
        ),
        "",
        (
            "`models/"
            "final_epigenomic_model_config.json`"
        ),
        "",
        "## Final workflow",
        "",
        "- DNA methylation preparation and QC",
        "- Patient matching",
        "- Multi-omics integration",
        "- Feature selection",
        "- Classification modelling",
        "- Survival modelling",
        "- Feature importance",
        "- Leakage-aware cross-validation",
        "- Threshold optimisation",
        "- Probability calibration",
        "- Final calibrated model",
        "",
        "## Important limitation",
        "",
        (
            "The model is a development-stage research model. "
            "Internal cross-validation does not replace "
            "independent external validation."
        ),
        "",
        (
            "The saved probabilities and classification "
            "threshold should not be used for clinical "
            "decision-making."
        ),
        "",
        "## Reproducibility",
        "",
        (
            "`manifest.csv` contains SHA256 checksums "
            "for packaged files."
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
    """Update final report and rebuild package."""

    print(
        "Updating final Month 2 report..."
    )

    updated_report = (
        create_updated_final_report()
    )

    FINAL_REPORT_FILE.write_text(
        updated_report,
        encoding="utf-8",
    )

    print(
        f"Updated report saved to:\n"
        f"{FINAL_REPORT_FILE}"
    )

    print(
        "\nRebuilding Month 2 package..."
    )

    prepare_package()

    reports_count = copy_all_files(
        source_dir=REPORTS_DIR,
        destination_dir=(
            PACKAGE_REPORTS_DIR
        ),
        extensions=(
            ".md",
            ".txt",
        ),
    )

    tables_count = copy_all_files(
        source_dir=TABLES_DIR,
        destination_dir=(
            PACKAGE_TABLES_DIR
        ),
        extensions=(
            ".csv",
            ".tsv",
        ),
    )

    figures_count = copy_all_files(
        source_dir=FIGURES_DIR,
        destination_dir=(
            PACKAGE_FIGURES_DIR
        ),
        extensions=(
            ".png",
            ".jpg",
            ".jpeg",
            ".svg",
            ".pdf",
        ),
    )

    models_count = copy_all_files(
        source_dir=MODELS_DIR,
        destination_dir=(
            PACKAGE_MODELS_DIR
        ),
        extensions=(
            ".pkl",
            ".joblib",
            ".json",
        ),
    )

    # =================================================
    # README
    # =================================================

    readme = create_package_readme(
        reports_count=reports_count,
        tables_count=tables_count,
        figures_count=figures_count,
        models_count=models_count,
    )

    PACKAGE_README_FILE.write_text(
        readme,
        encoding="utf-8",
    )

    # =================================================
    # MANIFEST
    # =================================================

    manifest = create_manifest()

    manifest.to_csv(
        PACKAGE_MANIFEST_FILE,
        index=False,
    )

    # =================================================
    # METADATA
    # =================================================

    metadata = {
        "project":
            "MI-Agata-MultiOmicsFlow",
        "package":
            "Month 2 final epigenomics package",
        "reports":
            reports_count,
        "tables":
            tables_count,
        "figures":
            figures_count,
        "models_and_configs":
            models_count,
        "final_model":
            (
                "models/"
                "final_epigenomic_calibrated_model.pkl"
            ),
        "calibration":
            "isotonic",
        "external_validation":
            False,
        "clinical_use":
            False,
    }

    if FINAL_MODEL_CONFIG_FILE.exists():

        try:
            final_config = json.loads(
                FINAL_MODEL_CONFIG_FILE.read_text(
                    encoding="utf-8"
                )
            )

            metadata[
                "final_model_configuration"
            ] = final_config

        except json.JSONDecodeError:

            metadata[
                "final_model_configuration"
            ] = (
                "Could not parse configuration JSON."
            )

    PACKAGE_METADATA_FILE.write_text(
        json.dumps(
            metadata,
            indent=2,
        ),
        encoding="utf-8",
    )

    # =================================================
    # DONE
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "FINAL PACKAGE COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        f"\nReports packaged: "
        f"{reports_count}"
    )

    print(
        f"Tables packaged: "
        f"{tables_count}"
    )

    print(
        f"Figures packaged: "
        f"{figures_count}"
    )

    print(
        f"Models/configs packaged: "
        f"{models_count}"
    )

    print(
        "\nFinal package saved to:\n"
        f"{PACKAGE_DIR}"
    )

    print(
        "\nFinal report:\n"
        f"{FINAL_REPORT_FILE}"
    )

    print(
        "\nREADME:\n"
        f"{PACKAGE_README_FILE}"
    )

    print(
        "\nManifest:\n"
        f"{PACKAGE_MANIFEST_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()