from pathlib import Path

import pandas as pd


# =====================================================
# PATHS
# =====================================================

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

REPORTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

RESULTS_SUMMARY_FILE = (
    REPORTS_DIR
    / "month2_results_summary.md"
)

KEY_METRICS_FILE = (
    TABLES_DIR
    / "month2_key_results.csv"
)


# =====================================================
# INPUT FILES
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


# =====================================================
# HELPERS
# =====================================================

def load_csv(
    file_path: Path,
) -> pd.DataFrame:
    """Load required CSV file."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def get_model_row(
    df: pd.DataFrame,
    model_name: str,
) -> pd.Series:
    """Return one model row."""

    rows = df[
        df["model"] == model_name
    ]

    if rows.empty:
        raise ValueError(
            f"Model not found in validation table: "
            f"{model_name}"
        )

    return rows.iloc[0]


def get_strategy_row(
    df: pd.DataFrame,
    strategy: str,
) -> pd.Series:
    """Return one threshold strategy row."""

    rows = df[
        df["strategy"] == strategy
    ]

    if rows.empty:
        raise ValueError(
            f"Strategy not found: {strategy}"
        )

    return rows.iloc[0]


def get_calibration_row(
    df: pd.DataFrame,
    method: str,
) -> pd.Series:
    """Return one calibration method row."""

    rows = df[
        df["method"] == method
    ]

    if rows.empty:
        raise ValueError(
            f"Calibration method not found: "
            f"{method}"
        )

    return rows.iloc[0]


def metric_value(
    df: pd.DataFrame,
    metric: str,
) -> float:
    """Read a metric from metric/value table."""

    rows = df[
        df["metric"] == metric
    ]

    if rows.empty:
        raise ValueError(
            f"Metric not found: {metric}"
        )

    return float(
        rows.iloc[0]["value"]
    )


def dataframe_to_markdown(
    df: pd.DataFrame,
) -> str:
    """Convert DataFrame to Markdown."""

    try:
        return df.to_markdown(
            index=False
        )

    except ImportError:
        return (
            "```text\n"
            + df.to_string(
                index=False
            )
            + "\n```"
        )


# =====================================================
# BUILD RESULTS
# =====================================================

def build_key_results() -> pd.DataFrame:
    """Build compact table of key Month 2 results."""

    validation = load_csv(
        VALIDATION_SUMMARY_FILE
    )

    threshold = load_csv(
        THRESHOLD_SUMMARY_FILE
    )

    calibration = load_csv(
        CALIBRATION_SUMMARY_FILE
    )

    final_metrics = load_csv(
        FINAL_METRICS_FILE
    )

    clinical = get_model_row(
        validation,
        "clinical_genomic",
    )

    epigenomic = get_model_row(
        validation,
        "clinical_genomic_epigenomic",
    )

    default_threshold = get_strategy_row(
        threshold,
        "default_0.50",
    )

    optimized_threshold = get_strategy_row(
        threshold,
        "optimized",
    )

    uncalibrated = get_calibration_row(
        calibration,
        "uncalibrated",
    )

    sigmoid = get_calibration_row(
        calibration,
        "sigmoid",
    )

    isotonic = get_calibration_row(
        calibration,
        "isotonic",
    )

    final_threshold = metric_value(
        final_metrics,
        "threshold",
    )

    final_roc_auc = metric_value(
        final_metrics,
        "roc_auc",
    )

    final_pr_auc = metric_value(
        final_metrics,
        "pr_auc",
    )

    final_brier = metric_value(
        final_metrics,
        "brier_score",
    )

    final_balanced_accuracy = metric_value(
        final_metrics,
        "balanced_accuracy",
    )

    final_precision = metric_value(
        final_metrics,
        "precision",
    )

    final_recall = metric_value(
        final_metrics,
        "recall",
    )

    final_f1 = metric_value(
        final_metrics,
        "f1_score",
    )

    rows = [
        {
            "result":
                "Clinical-genomic ROC AUC",
            "value":
                clinical["roc_auc_mean"],
            "std":
                clinical["roc_auc_std"],
        },
        {
            "result":
                "Clinical-genomic-epigenomic ROC AUC",
            "value":
                epigenomic["roc_auc_mean"],
            "std":
                epigenomic["roc_auc_std"],
        },
        {
            "result":
                "ROC AUC improvement after methylation",
            "value":
                (
                    epigenomic["roc_auc_mean"]
                    - clinical["roc_auc_mean"]
                ),
            "std":
                None,
        },
        {
            "result":
                "Default threshold recall",
            "value":
                default_threshold[
                    "recall_mean"
                ],
            "std":
                default_threshold[
                    "recall_std"
                ],
        },
        {
            "result":
                "Optimized threshold recall",
            "value":
                optimized_threshold[
                    "recall_mean"
                ],
            "std":
                optimized_threshold[
                    "recall_std"
                ],
        },
        {
            "result":
                "Default threshold balanced accuracy",
            "value":
                default_threshold[
                    "balanced_accuracy_mean"
                ],
            "std":
                default_threshold[
                    "balanced_accuracy_std"
                ],
        },
        {
            "result":
                "Optimized threshold balanced accuracy",
            "value":
                optimized_threshold[
                    "balanced_accuracy_mean"
                ],
            "std":
                optimized_threshold[
                    "balanced_accuracy_std"
                ],
        },
        {
            "result":
                "Optimized uncalibrated threshold",
            "value":
                optimized_threshold[
                    "threshold_mean"
                ],
            "std":
                optimized_threshold[
                    "threshold_std"
                ],
        },
        {
            "result":
                "Uncalibrated Brier score",
            "value":
                uncalibrated[
                    "brier_score_mean"
                ],
            "std":
                uncalibrated[
                    "brier_score_std"
                ],
        },
        {
            "result":
                "Sigmoid Brier score",
            "value":
                sigmoid[
                    "brier_score_mean"
                ],
            "std":
                sigmoid[
                    "brier_score_std"
                ],
        },
        {
            "result":
                "Isotonic Brier score",
            "value":
                isotonic[
                    "brier_score_mean"
                ],
            "std":
                isotonic[
                    "brier_score_std"
                ],
        },
        {
            "result":
                "Final calibrated ROC AUC",
            "value":
                final_roc_auc,
            "std":
                None,
        },
        {
            "result":
                "Final calibrated PR AUC",
            "value":
                final_pr_auc,
            "std":
                None,
        },
        {
            "result":
                "Final calibrated Brier score",
            "value":
                final_brier,
            "std":
                None,
        },
        {
            "result":
                "Final calibrated threshold",
            "value":
                final_threshold,
            "std":
                None,
        },
        {
            "result":
                "Final balanced accuracy",
            "value":
                final_balanced_accuracy,
            "std":
                None,
        },
        {
            "result":
                "Final precision",
            "value":
                final_precision,
            "std":
                None,
        },
        {
            "result":
                "Final recall",
            "value":
                final_recall,
            "std":
                None,
        },
        {
            "result":
                "Final F1-score",
            "value":
                final_f1,
            "std":
                None,
        },
    ]

    return pd.DataFrame(
        rows
    )


# =====================================================
# PUBLICATION SUMMARY
# =====================================================

def create_results_summary(
    key_results: pd.DataFrame,
) -> str:
    """Create publication-style Results summary."""

    validation = load_csv(
        VALIDATION_SUMMARY_FILE
    )

    threshold = load_csv(
        THRESHOLD_SUMMARY_FILE
    )

    calibration = load_csv(
        CALIBRATION_SUMMARY_FILE
    )

    final_metrics = load_csv(
        FINAL_METRICS_FILE
    )

    clinical = get_model_row(
        validation,
        "clinical_genomic",
    )

    epigenomic = get_model_row(
        validation,
        "clinical_genomic_epigenomic",
    )

    default_threshold = get_strategy_row(
        threshold,
        "default_0.50",
    )

    optimized_threshold = get_strategy_row(
        threshold,
        "optimized",
    )

    isotonic = get_calibration_row(
        calibration,
        "isotonic",
    )

    uncalibrated = get_calibration_row(
        calibration,
        "uncalibrated",
    )

    final_threshold = metric_value(
        final_metrics,
        "threshold",
    )

    final_roc_auc = metric_value(
        final_metrics,
        "roc_auc",
    )

    final_pr_auc = metric_value(
        final_metrics,
        "pr_auc",
    )

    final_brier = metric_value(
        final_metrics,
        "brier_score",
    )

    final_balanced_accuracy = metric_value(
        final_metrics,
        "balanced_accuracy",
    )

    final_precision = metric_value(
        final_metrics,
        "precision",
    )

    final_recall = metric_value(
        final_metrics,
        "recall",
    )

    final_f1 = metric_value(
        final_metrics,
        "f1_score",
    )

    auc_improvement = (
        epigenomic["roc_auc_mean"]
        - clinical["roc_auc_mean"]
    )

    report_lines = [
        "# Month 2 Results Summary",
        "",
        "## Results",
        "",
        (
            "Integration of DNA methylation substantially "
            "improved predictive discrimination compared with "
            "the clinical-genomic baseline."
        ),
        "",
        (
            "In leakage-aware five-fold cross-validation, "
            "the clinical-genomic model achieved a mean ROC AUC "
            f"of {clinical['roc_auc_mean']:.3f} "
            f"± {clinical['roc_auc_std']:.3f}, whereas the "
            "clinical-genomic-epigenomic model achieved "
            f"{epigenomic['roc_auc_mean']:.3f} "
            f"± {epigenomic['roc_auc_std']:.3f}."
        ),
        "",
        (
            "This corresponded to an absolute ROC AUC increase "
            f"of {auc_improvement:.3f} after inclusion of "
            "methylation features."
        ),
        "",
        "## Classification threshold",
        "",
        (
            "At the conventional probability threshold of 0.50, "
            "the epigenomic classifier showed high overall "
            "accuracy but low sensitivity for the positive class."
        ),
        "",
        (
            f"Mean recall was "
            f"{default_threshold['recall_mean']:.3f}, "
            "with balanced accuracy of "
            f"{default_threshold['balanced_accuracy_mean']:.3f}."
        ),
        "",
        (
            "Training-derived threshold optimisation improved "
            "minority-class detection."
        ),
        "",
        (
            f"The mean optimized uncalibrated threshold was "
            f"{optimized_threshold['threshold_mean']:.3f} "
            f"± {optimized_threshold['threshold_std']:.3f}."
        ),
        "",
        (
            f"At the optimized threshold, recall increased to "
            f"{optimized_threshold['recall_mean']:.3f}, "
            "balanced accuracy increased to "
            f"{optimized_threshold['balanced_accuracy_mean']:.3f}, "
            "and F1-score increased to "
            f"{optimized_threshold['f1_score_mean']:.3f}."
        ),
        "",
        "## Probability calibration",
        "",
        (
            "Probability calibration improved the quality of "
            "predicted risk estimates."
        ),
        "",
        (
            f"The uncalibrated model had a mean Brier score of "
            f"{uncalibrated['brier_score_mean']:.4f}, whereas "
            "isotonic calibration reduced the mean Brier score "
            f"to {isotonic['brier_score_mean']:.4f}."
        ),
        "",
        (
            f"Isotonic calibration achieved a mean ROC AUC of "
            f"{isotonic['roc_auc_mean']:.3f} and mean PR AUC of "
            f"{isotonic['pr_auc_mean']:.3f}."
        ),
        "",
        "## Final development model",
        "",
        (
            "The final development-stage classifier consisted "
            "of a Random Forest model with isotonic probability "
            "calibration."
        ),
        "",
        (
            f"The final calibrated operating threshold was "
            f"{final_threshold:.3f}."
        ),
        "",
        (
            f"Leakage-aware out-of-fold performance of the "
            f"final model included ROC AUC {final_roc_auc:.3f}, "
            f"PR AUC {final_pr_auc:.3f}, "
            f"Brier score {final_brier:.4f}, "
            f"balanced accuracy {final_balanced_accuracy:.3f}, "
            f"precision {final_precision:.3f}, "
            f"recall {final_recall:.3f}, and "
            f"F1-score {final_f1:.3f}."
        ),
        "",
        "## Key results table",
        "",
        dataframe_to_markdown(
            key_results.round(4)
        ),
        "",
        "## Recommended figures",
        "",
        "1. `model_validation_roc_auc.png`",
        "2. `threshold_validation_roc_curve.png`",
        "3. `threshold_validation_precision_recall.png`",
        "4. `confusion_matrix_optimized_threshold.png`",
        "5. `model_calibration_curve.png`",
        "6. `cox_risk_groups_survival.png`",
        "7. `cox_top_hazard_ratios.png`",
        "8. `combined_feature_importance.png`",
        "",
        "## Interpretation",
        "",
        (
            "The results support the hypothesis that DNA "
            "methylation contains predictive information beyond "
            "the clinical-genomic baseline."
        ),
        "",
        (
            "The decline from the original single-split ROC AUC "
            "of 1.0 to approximately 0.86 under leakage-aware "
            "cross-validation demonstrates the importance of "
            "rigorous validation."
        ),
        "",
        (
            "Threshold optimisation substantially improved "
            "positive-class sensitivity, while probability "
            "calibration improved agreement between predicted "
            "probabilities and observed outcomes."
        ),
        "",
        "## Limitations",
        "",
        (
            "All reported model-performance estimates are based "
            "on internal validation within the available cohort."
        ),
        "",
        (
            "External validation on an independent dataset is "
            "required before assessing transportability or "
            "clinical utility."
        ),
        "",
        (
            "The optimized classification threshold is a "
            "development-stage operating point and should not "
            "be interpreted as a validated clinical cutoff."
        ),
        "",
    ]

    return "\n".join(
        report_lines
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Generate publication-style Month 2 results summary."""

    print(
        "Generating Month 2 results summary..."
    )

    key_results = build_key_results()

    key_results.to_csv(
        KEY_METRICS_FILE,
        index=False,
    )

    report = create_results_summary(
        key_results
    )

    RESULTS_SUMMARY_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print(
        "\n========================================"
    )

    print(
        "RESULTS SUMMARY COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        "\nKey results:"
    )

    print(
        key_results.to_string(
            index=False
        )
    )

    print(
        "\nResults summary saved to:\n"
        f"{RESULTS_SUMMARY_FILE}"
    )

    print(
        "\nKey results table saved to:\n"
        f"{KEY_METRICS_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()