from pathlib import Path
import json
import math

import pandas as pd


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUTS_DIR = PROJECT_ROOT / "outputs"

TABLES_DIR = OUTPUTS_DIR / "tables"
FIGURES_DIR = OUTPUTS_DIR / "figures"
REPORTS_DIR = OUTPUTS_DIR / "reports"
MODELS_DIR = OUTPUTS_DIR / "models"

PACKAGE_DIR = OUTPUTS_DIR / "month2_package"

QC_TABLE_FILE = (
    TABLES_DIR
    / "month2_final_qc.csv"
)

QC_REPORT_FILE = (
    REPORTS_DIR
    / "month2_final_qc_report.md"
)


# =====================================================
# REQUIRED FILES
# =====================================================

REQUIRED_REPORTS = [
    "month2_final_report.md",
    "model_validation_report.md",
    "threshold_optimization_report.md",
    "model_calibration_report.md",
    "final_epigenomic_model_report.md",
    "month2_results_summary.md",
    "month2_manuscript_results.md",
    "month2_manuscript_discussion.md",
    "month2_manuscript_conclusions.md",
    "month2_manuscript_sections.md",
]

REQUIRED_TABLES = [
    "model_validation_summary.csv",
    "model_validation_folds.csv",
    "model_leakage_audit.csv",
    "threshold_optimization_summary.csv",
    "optimized_thresholds.csv",
    "calibration_validation_summary.csv",
    "calibration_validation_folds.csv",
    "final_epigenomic_model_metrics.csv",
    "final_epigenomic_model_features.csv",
    "final_calibrated_threshold.csv",
    "final_calibrated_oof_predictions.csv",
    "month2_key_results.csv",
]

REQUIRED_FIGURES = [
    "model_validation_roc_auc.png",
    "threshold_validation_roc_curve.png",
    "threshold_validation_precision_recall.png",
    "confusion_matrix_threshold_050.png",
    "confusion_matrix_optimized_threshold.png",
    "model_calibration_curve.png",
    "calibrated_probability_distribution.png",
    "cox_risk_groups_survival.png",
    "cox_top_hazard_ratios.png",
    "combined_feature_importance.png",
]

REQUIRED_MODELS = [
    "final_epigenomic_calibrated_model.pkl",
    "final_epigenomic_model_config.json",
]


# =====================================================
# HELPERS
# =====================================================

def add_check(
    checks: list[dict],
    category: str,
    check: str,
    status: str,
    details: str,
) -> None:
    """Add one QC result."""

    checks.append(
        {
            "category": category,
            "check": check,
            "status": status,
            "details": details,
        }
    )


def check_required_files(
    checks: list[dict],
    directory: Path,
    filenames: list[str],
    category: str,
) -> None:
    """Check whether required files exist."""

    for filename in filenames:

        file_path = directory / filename

        if file_path.exists():

            size = file_path.stat().st_size

            if size > 0:

                add_check(
                    checks,
                    category,
                    filename,
                    "PASS",
                    f"Exists; size={size:,} bytes",
                )

            else:

                add_check(
                    checks,
                    category,
                    filename,
                    "FAIL",
                    "File exists but is empty",
                )

        else:

            add_check(
                checks,
                category,
                filename,
                "FAIL",
                f"Missing: {file_path}",
            )


def load_csv_safe(
    file_path: Path,
) -> pd.DataFrame | None:
    """Load CSV safely."""

    if not file_path.exists():
        return None

    try:
        return pd.read_csv(
            file_path,
            low_memory=False,
        )

    except Exception:
        return None


def metric_value(
    df: pd.DataFrame,
    metric: str,
) -> float | None:
    """Extract metric from metric/value table."""

    if (
        "metric" not in df.columns
        or "value" not in df.columns
    ):
        return None

    rows = df[
        df["metric"] == metric
    ]

    if rows.empty:
        return None

    try:
        return float(
            rows.iloc[0]["value"]
        )

    except (TypeError, ValueError):
        return None


def value_is_finite(
    value,
) -> bool:
    """Check whether numeric value is finite."""

    try:
        return math.isfinite(
            float(value)
        )

    except (TypeError, ValueError):
        return False


# =====================================================
# FINAL METRICS QC
# =====================================================

def check_final_metrics(
    checks: list[dict],
) -> None:
    """Validate final metric ranges."""

    file_path = (
        TABLES_DIR
        / "final_epigenomic_model_metrics.csv"
    )

    df = load_csv_safe(
        file_path
    )

    if df is None:

        add_check(
            checks,
            "metrics",
            "Final metrics readable",
            "FAIL",
            "Could not load final metrics CSV",
        )

        return

    required_metrics = {
        "roc_auc": (0.0, 1.0),
        "pr_auc": (0.0, 1.0),
        "brier_score": (0.0, 1.0),
        "threshold": (0.0, 1.0),
        "accuracy": (0.0, 1.0),
        "balanced_accuracy": (0.0, 1.0),
        "precision": (0.0, 1.0),
        "recall": (0.0, 1.0),
        "f1_score": (0.0, 1.0),
    }

    for metric, (
        minimum,
        maximum,
    ) in required_metrics.items():

        value = metric_value(
            df,
            metric,
        )

        if value is None:

            add_check(
                checks,
                "metrics",
                metric,
                "FAIL",
                "Metric missing",
            )

            continue

        if not value_is_finite(
            value
        ):

            add_check(
                checks,
                "metrics",
                metric,
                "FAIL",
                f"Non-finite value: {value}",
            )

            continue

        if not (
            minimum
            <= value
            <= maximum
        ):

            add_check(
                checks,
                "metrics",
                metric,
                "FAIL",
                (
                    f"Out of expected range "
                    f"[{minimum}, {maximum}]: "
                    f"{value}"
                ),
            )

        else:

            add_check(
                checks,
                "metrics",
                metric,
                "PASS",
                f"{value:.6f}",
            )


# =====================================================
# VALIDATION QC
# =====================================================

def check_validation_summary(
    checks: list[dict],
) -> None:
    """Check Step 21 validation output."""

    file_path = (
        TABLES_DIR
        / "model_validation_summary.csv"
    )

    df = load_csv_safe(
        file_path
    )

    if df is None:

        add_check(
            checks,
            "validation",
            "Validation summary readable",
            "FAIL",
            "Could not load model validation summary",
        )

        return

    expected_models = {
        "clinical_genomic",
        "clinical_genomic_epigenomic",
    }

    if "model" not in df.columns:

        add_check(
            checks,
            "validation",
            "Model column",
            "FAIL",
            "Column 'model' missing",
        )

        return

    actual_models = set(
        df["model"].astype(str)
    )

    missing = (
        expected_models
        - actual_models
    )

    if missing:

        add_check(
            checks,
            "validation",
            "Expected validation models",
            "FAIL",
            f"Missing models: {sorted(missing)}",
        )

    else:

        add_check(
            checks,
            "validation",
            "Expected validation models",
            "PASS",
            "Both baseline and epigenomic models present",
        )

    if "roc_auc_mean" in df.columns:

        for _, row in df.iterrows():

            value = row[
                "roc_auc_mean"
            ]

            if value_is_finite(value) and (
                0 <= float(value) <= 1
            ):

                add_check(
                    checks,
                    "validation",
                    (
                        f"ROC AUC: "
                        f"{row['model']}"
                    ),
                    "PASS",
                    f"{float(value):.6f}",
                )

            else:

                add_check(
                    checks,
                    "validation",
                    (
                        f"ROC AUC: "
                        f"{row['model']}"
                    ),
                    "FAIL",
                    f"Invalid value: {value}",
                )


# =====================================================
# THRESHOLD QC
# =====================================================

def check_threshold_results(
    checks: list[dict],
) -> None:
    """Check threshold optimization."""

    file_path = (
        TABLES_DIR
        / "threshold_optimization_summary.csv"
    )

    df = load_csv_safe(
        file_path
    )

    if df is None:

        add_check(
            checks,
            "threshold",
            "Threshold summary readable",
            "FAIL",
            "Could not load threshold summary",
        )

        return

    if "strategy" not in df.columns:

        add_check(
            checks,
            "threshold",
            "Strategy column",
            "FAIL",
            "Column 'strategy' missing",
        )

        return

    strategies = set(
        df["strategy"].astype(str)
    )

    expected = {
        "default_0.50",
        "optimized",
    }

    if expected.issubset(
        strategies
    ):

        add_check(
            checks,
            "threshold",
            "Threshold strategies",
            "PASS",
            "Default and optimized strategies present",
        )

    else:

        add_check(
            checks,
            "threshold",
            "Threshold strategies",
            "FAIL",
            f"Found strategies: {sorted(strategies)}",
        )

    optimized = df[
        df["strategy"]
        == "optimized"
    ]

    if not optimized.empty:

        value = optimized.iloc[0].get(
            "threshold_mean"
        )

        if value_is_finite(
            value
        ) and (
            0 < float(value) < 1
        ):

            add_check(
                checks,
                "threshold",
                "Optimized raw threshold",
                "PASS",
                f"{float(value):.6f}",
            )

        else:

            add_check(
                checks,
                "threshold",
                "Optimized raw threshold",
                "FAIL",
                f"Invalid threshold: {value}",
            )


# =====================================================
# CALIBRATION QC
# =====================================================

def check_calibration_results(
    checks: list[dict],
) -> None:
    """Check probability calibration."""

    file_path = (
        TABLES_DIR
        / "calibration_validation_summary.csv"
    )

    df = load_csv_safe(
        file_path
    )

    if df is None:

        add_check(
            checks,
            "calibration",
            "Calibration summary readable",
            "FAIL",
            "Could not load calibration summary",
        )

        return

    if "method" not in df.columns:

        add_check(
            checks,
            "calibration",
            "Calibration method column",
            "FAIL",
            "Column 'method' missing",
        )

        return

    expected = {
        "uncalibrated",
        "sigmoid",
        "isotonic",
    }

    actual = set(
        df["method"].astype(str)
    )

    if expected.issubset(
        actual
    ):

        add_check(
            checks,
            "calibration",
            "Calibration methods",
            "PASS",
            "All three methods present",
        )

    else:

        add_check(
            checks,
            "calibration",
            "Calibration methods",
            "FAIL",
            f"Found: {sorted(actual)}",
        )

    if "brier_score_mean" in df.columns:

        valid = df.dropna(
            subset=[
                "brier_score_mean"
            ]
        )

        if not valid.empty:

            best = valid.sort_values(
                "brier_score_mean"
            ).iloc[0]

            best_method = str(
                best["method"]
            )

            if best_method == "isotonic":

                add_check(
                    checks,
                    "calibration",
                    "Best calibration by Brier score",
                    "PASS",
                    (
                        f"isotonic; "
                        f"Brier="
                        f"{best['brier_score_mean']:.6f}"
                    ),
                )

            else:

                add_check(
                    checks,
                    "calibration",
                    "Best calibration by Brier score",
                    "WARNING",
                    (
                        f"Best method is "
                        f"{best_method}, "
                        "not isotonic"
                    ),
                )


# =====================================================
# CONFIG CONSISTENCY
# =====================================================

def check_model_config(
    checks: list[dict],
) -> None:
    """Compare config with final metrics."""

    config_file = (
        MODELS_DIR
        / "final_epigenomic_model_config.json"
    )

    metrics_file = (
        TABLES_DIR
        / "final_epigenomic_model_metrics.csv"
    )

    if not config_file.exists():

        add_check(
            checks,
            "config",
            "Final config exists",
            "FAIL",
            "Configuration JSON missing",
        )

        return

    try:

        config = json.loads(
            config_file.read_text(
                encoding="utf-8"
            )
        )

    except Exception as error:

        add_check(
            checks,
            "config",
            "Final config readable",
            "FAIL",
            f"{type(error).__name__}: {error}",
        )

        return

    add_check(
        checks,
        "config",
        "Final config readable",
        "PASS",
        "Valid JSON",
    )

    calibration = config.get(
        "calibration"
    )

    if calibration == "isotonic":

        add_check(
            checks,
            "config",
            "Calibration configuration",
            "PASS",
            "isotonic",
        )

    else:

        add_check(
            checks,
            "config",
            "Calibration configuration",
            "WARNING",
            f"Configured as: {calibration}",
        )

    metrics = load_csv_safe(
        metrics_file
    )

    if metrics is None:
        return

    metric_threshold = metric_value(
        metrics,
        "threshold",
    )

    config_threshold = config.get(
        "development_threshold"
    )

    if (
        metric_threshold is not None
        and config_threshold is not None
    ):

        difference = abs(
            float(metric_threshold)
            - float(config_threshold)
        )

        if difference < 1e-9:

            add_check(
                checks,
                "config",
                "Threshold consistency",
                "PASS",
                (
                    f"Config and metrics agree: "
                    f"{metric_threshold:.6f}"
                ),
            )

        else:

            add_check(
                checks,
                "config",
                "Threshold consistency",
                "FAIL",
                (
                    f"Metrics={metric_threshold}; "
                    f"config={config_threshold}"
                ),
            )


# =====================================================
# OOF PREDICTIONS
# =====================================================

def check_oof_predictions(
    checks: list[dict],
) -> None:
    """Validate final OOF prediction table."""

    file_path = (
        TABLES_DIR
        / "final_calibrated_oof_predictions.csv"
    )

    df = load_csv_safe(
        file_path
    )

    if df is None:

        add_check(
            checks,
            "OOF",
            "OOF predictions readable",
            "FAIL",
            "Could not load final OOF predictions",
        )

        return

    required_columns = {
        "fold",
        "patient_id",
        "true_outcome",
        "calibrated_probability",
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:

        add_check(
            checks,
            "OOF",
            "Required OOF columns",
            "FAIL",
            f"Missing: {sorted(missing)}",
        )

        return

    add_check(
        checks,
        "OOF",
        "Required OOF columns",
        "PASS",
        "All required columns present",
    )

    duplicate_patients = int(
        df["patient_id"]
        .duplicated()
        .sum()
    )

    if duplicate_patients == 0:

        add_check(
            checks,
            "OOF",
            "Unique OOF patients",
            "PASS",
            f"{len(df):,} unique patients",
        )

    else:

        add_check(
            checks,
            "OOF",
            "Unique OOF patients",
            "FAIL",
            (
                f"{duplicate_patients:,} "
                "duplicated patient IDs"
            ),
        )

    probabilities = pd.to_numeric(
        df["calibrated_probability"],
        errors="coerce",
    )

    missing_probabilities = int(
        probabilities.isna().sum()
    )

    out_of_range = int(
        (
            (probabilities < 0)
            | (probabilities > 1)
        )
        .sum()
    )

    if (
        missing_probabilities == 0
        and out_of_range == 0
    ):

        add_check(
            checks,
            "OOF",
            "Probability validity",
            "PASS",
            "All probabilities are finite and in [0,1]",
        )

    else:

        add_check(
            checks,
            "OOF",
            "Probability validity",
            "FAIL",
            (
                f"NaN={missing_probabilities}; "
                f"out_of_range={out_of_range}"
            ),
        )


# =====================================================
# PACKAGE QC
# =====================================================

def check_package(
    checks: list[dict],
) -> None:
    """Check final package completeness."""

    if not PACKAGE_DIR.exists():

        add_check(
            checks,
            "package",
            "Month 2 package",
            "FAIL",
            f"Missing: {PACKAGE_DIR}",
        )

        return

    add_check(
        checks,
        "package",
        "Month 2 package",
        "PASS",
        str(PACKAGE_DIR),
    )

    package_required = [
        PACKAGE_DIR / "README.md",
        PACKAGE_DIR / "manifest.csv",
        PACKAGE_DIR / "package_metadata.json",
        (
            PACKAGE_DIR
            / "models"
            / "final_epigenomic_calibrated_model.pkl"
        ),
        (
            PACKAGE_DIR
            / "models"
            / "final_epigenomic_model_config.json"
        ),
        (
            PACKAGE_DIR
            / "reports"
            / "month2_final_report.md"
        ),
    ]

    for file_path in package_required:

        if (
            file_path.exists()
            and file_path.stat().st_size > 0
        ):

            add_check(
                checks,
                "package",
                (
                    "Package file: "
                    f"{file_path.name}"
                ),
                "PASS",
                str(
                    file_path.relative_to(
                        PACKAGE_DIR
                    )
                ),
            )

        else:

            add_check(
                checks,
                "package",
                (
                    "Package file: "
                    f"{file_path.name}"
                ),
                "FAIL",
                f"Missing or empty: {file_path}",
            )


# =====================================================
# MANUSCRIPT QC
# =====================================================

def check_manuscript_files(
    checks: list[dict],
) -> None:
    """Check manuscript sections are not empty."""

    manuscript_files = [
        (
            REPORTS_DIR
            / "month2_manuscript_results.md"
        ),
        (
            REPORTS_DIR
            / "month2_manuscript_discussion.md"
        ),
        (
            REPORTS_DIR
            / "month2_manuscript_conclusions.md"
        ),
        (
            REPORTS_DIR
            / "month2_manuscript_sections.md"
        ),
    ]

    for file_path in manuscript_files:

        if not file_path.exists():

            add_check(
                checks,
                "manuscript",
                file_path.name,
                "FAIL",
                "File missing",
            )

            continue

        text = file_path.read_text(
            encoding="utf-8"
        ).strip()

        if len(text) >= 100:

            add_check(
                checks,
                "manuscript",
                file_path.name,
                "PASS",
                (
                    f"{len(text):,} characters"
                ),
            )

        else:

            add_check(
                checks,
                "manuscript",
                file_path.name,
                "WARNING",
                (
                    f"Very short file: "
                    f"{len(text):,} characters"
                ),
            )


# =====================================================
# REPORT
# =====================================================

def create_qc_report(
    qc: pd.DataFrame,
    overall_status: str,
) -> str:
    """Generate Markdown QC report."""

    status_counts = (
        qc["status"]
        .value_counts()
        .to_dict()
    )

    failures = qc[
        qc["status"] == "FAIL"
    ]

    warnings = qc[
        qc["status"] == "WARNING"
    ]

    try:

        table_text = qc.to_markdown(
            index=False
        )

    except ImportError:

        table_text = (
            "```text\n"
            + qc.to_string(
                index=False
            )
            + "\n```"
        )

    lines = [
        "# Month 2 Final Quality-Control Report",
        "",
        "## Overall status",
        "",
        f"**{overall_status}**",
        "",
        "## Summary",
        "",
        (
            f"- PASS: "
            f"{status_counts.get('PASS', 0)}"
        ),
        (
            f"- WARNING: "
            f"{status_counts.get('WARNING', 0)}"
        ),
        (
            f"- FAIL: "
            f"{status_counts.get('FAIL', 0)}"
        ),
        "",
        "## Interpretation",
        "",
    ]

    if overall_status == "PASS":

        lines.extend(
            [
                (
                    "All critical Month 2 quality-control "
                    "checks passed."
                ),
                "",
                (
                    "The analytical pipeline and final "
                    "development package are internally "
                    "consistent according to the automated "
                    "checks implemented in this audit."
                ),
            ]
        )

    elif overall_status == "WARNING":

        lines.extend(
            [
                (
                    "No critical failure was detected, but "
                    "one or more non-critical warnings "
                    "require review."
                ),
            ]
        )

    else:

        lines.extend(
            [
                (
                    "One or more critical quality-control "
                    "checks failed. The affected files or "
                    "results should be corrected before "
                    "the pipeline is considered final."
                ),
            ]
        )

    if not failures.empty:

        lines.extend(
            [
                "",
                "## Failed checks",
                "",
            ]
        )

        for _, row in failures.iterrows():

            lines.append(
                (
                    f"- **{row['category']} / "
                    f"{row['check']}**: "
                    f"{row['details']}"
                )
            )

    if not warnings.empty:

        lines.extend(
            [
                "",
                "## Warnings",
                "",
            ]
        )

        for _, row in warnings.iterrows():

            lines.append(
                (
                    f"- **{row['category']} / "
                    f"{row['check']}**: "
                    f"{row['details']}"
                )
            )

    lines.extend(
        [
            "",
            "## Complete QC table",
            "",
            table_text,
            "",
            "## Important limitation",
            "",
            (
                "Passing this automated QC confirms internal "
                "file and result consistency only. It does not "
                "constitute independent external validation "
                "of the predictive model."
            ),
            "",
        ]
    )

    return "\n".join(
        lines
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Run final Month 2 quality-control audit."""

    print(
        "Running final Month 2 QC audit..."
    )

    checks = []

    # ---------------------------------
    # FILES
    # ---------------------------------

    check_required_files(
        checks,
        REPORTS_DIR,
        REQUIRED_REPORTS,
        "required_report",
    )

    check_required_files(
        checks,
        TABLES_DIR,
        REQUIRED_TABLES,
        "required_table",
    )

    check_required_files(
        checks,
        FIGURES_DIR,
        REQUIRED_FIGURES,
        "required_figure",
    )

    check_required_files(
        checks,
        MODELS_DIR,
        REQUIRED_MODELS,
        "required_model",
    )

    # ---------------------------------
    # CONTENT / CONSISTENCY
    # ---------------------------------

    check_final_metrics(
        checks
    )

    check_validation_summary(
        checks
    )

    check_threshold_results(
        checks
    )

    check_calibration_results(
        checks
    )

    check_model_config(
        checks
    )

    check_oof_predictions(
        checks
    )

    check_package(
        checks
    )

    check_manuscript_files(
        checks
    )

    # ---------------------------------
    # RESULT
    # ---------------------------------

    qc = pd.DataFrame(
        checks
    )

    if (
        qc["status"] == "FAIL"
    ).any():

        overall_status = "FAIL"

    elif (
        qc["status"] == "WARNING"
    ).any():

        overall_status = "WARNING"

    else:

        overall_status = "PASS"

    # ---------------------------------
    # SAVE
    # ---------------------------------

    qc.to_csv(
        QC_TABLE_FILE,
        index=False,
    )

    report = create_qc_report(
        qc=qc,
        overall_status=overall_status,
    )

    QC_REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # ---------------------------------
    # TERMINAL
    # ---------------------------------

    print(
        "\n========================================"
    )

    print(
        "FINAL MONTH 2 QC"
    )

    print(
        "========================================"
    )

    print(
        f"\nFINAL QC STATUS: "
        f"{overall_status}"
    )

    print(
        "\nStatus counts:"
    )

    print(
        qc["status"]
        .value_counts()
        .to_string()
    )

    failures = qc[
        qc["status"] == "FAIL"
    ]

    warnings = qc[
        qc["status"] == "WARNING"
    ]

    if not failures.empty:

        print(
            "\nFAILED CHECKS:"
        )

        for _, row in failures.iterrows():

            print(
                f"- {row['category']} / "
                f"{row['check']}: "
                f"{row['details']}"
            )

    if not warnings.empty:

        print(
            "\nWARNINGS:"
        )

        for _, row in warnings.iterrows():

            print(
                f"- {row['category']} / "
                f"{row['check']}: "
                f"{row['details']}"
            )

    print(
        "\nQC table saved to:\n"
        f"{QC_TABLE_FILE}"
    )

    print(
        "\nQC report saved to:\n"
        f"{QC_REPORT_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()