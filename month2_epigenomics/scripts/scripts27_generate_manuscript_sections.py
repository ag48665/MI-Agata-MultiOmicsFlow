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


VALIDATION_FILE = (
    TABLES_DIR
    / "model_validation_summary.csv"
)

THRESHOLD_FILE = (
    TABLES_DIR
    / "threshold_optimization_summary.csv"
)

CALIBRATION_FILE = (
    TABLES_DIR
    / "calibration_validation_summary.csv"
)

FINAL_METRICS_FILE = (
    TABLES_DIR
    / "final_epigenomic_model_metrics.csv"
)

KEY_RESULTS_FILE = (
    TABLES_DIR
    / "month2_key_results.csv"
)


RESULTS_FILE = (
    REPORTS_DIR
    / "month2_manuscript_results.md"
)

DISCUSSION_FILE = (
    REPORTS_DIR
    / "month2_manuscript_discussion.md"
)

CONCLUSIONS_FILE = (
    REPORTS_DIR
    / "month2_manuscript_conclusions.md"
)

COMBINED_FILE = (
    REPORTS_DIR
    / "month2_manuscript_sections.md"
)


# =====================================================
# HELPERS
# =====================================================

def load_csv(
    file_path: Path,
) -> pd.DataFrame:
    """Load a required CSV file."""

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
    """Return one validation model row."""

    rows = df[
        df["model"] == model_name
    ]

    if rows.empty:
        raise ValueError(
            f"Model not found: {model_name}"
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
            f"Threshold strategy not found: "
            f"{strategy}"
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
    """Read one value from metric/value table."""

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


# =====================================================
# LOAD RESULTS
# =====================================================

def load_all_results() -> dict:
    """Load all result tables needed for manuscript text."""

    validation = load_csv(
        VALIDATION_FILE
    )

    threshold = load_csv(
        THRESHOLD_FILE
    )

    calibration = load_csv(
        CALIBRATION_FILE
    )

    final_metrics = load_csv(
        FINAL_METRICS_FILE
    )

    key_results = load_csv(
        KEY_RESULTS_FILE
    )

    clinical = get_model_row(
        validation,
        "clinical_genomic",
    )

    epigenomic = get_model_row(
        validation,
        "clinical_genomic_epigenomic",
    )

    default_threshold = (
        get_strategy_row(
            threshold,
            "default_0.50",
        )
    )

    optimized_threshold = (
        get_strategy_row(
            threshold,
            "optimized",
        )
    )

    uncalibrated = (
        get_calibration_row(
            calibration,
            "uncalibrated",
        )
    )

    sigmoid = (
        get_calibration_row(
            calibration,
            "sigmoid",
        )
    )

    isotonic = (
        get_calibration_row(
            calibration,
            "isotonic",
        )
    )

    final = {
        "roc_auc":
            metric_value(
                final_metrics,
                "roc_auc",
            ),
        "pr_auc":
            metric_value(
                final_metrics,
                "pr_auc",
            ),
        "brier_score":
            metric_value(
                final_metrics,
                "brier_score",
            ),
        "log_loss":
            metric_value(
                final_metrics,
                "log_loss",
            ),
        "threshold":
            metric_value(
                final_metrics,
                "threshold",
            ),
        "accuracy":
            metric_value(
                final_metrics,
                "accuracy",
            ),
        "balanced_accuracy":
            metric_value(
                final_metrics,
                "balanced_accuracy",
            ),
        "precision":
            metric_value(
                final_metrics,
                "precision",
            ),
        "recall":
            metric_value(
                final_metrics,
                "recall",
            ),
        "f1_score":
            metric_value(
                final_metrics,
                "f1_score",
            ),
    }

    return {
        "validation":
            validation,
        "threshold":
            threshold,
        "calibration":
            calibration,
        "final_metrics":
            final_metrics,
        "key_results":
            key_results,
        "clinical":
            clinical,
        "epigenomic":
            epigenomic,
        "default_threshold":
            default_threshold,
        "optimized_threshold":
            optimized_threshold,
        "uncalibrated":
            uncalibrated,
        "sigmoid":
            sigmoid,
        "isotonic":
            isotonic,
        "final":
            final,
    }


# =====================================================
# RESULTS SECTION
# =====================================================

def create_results_section(
    r: dict,
) -> str:
    """Create manuscript-style Results section."""

    clinical = r["clinical"]
    epigenomic = r["epigenomic"]

    default_threshold = (
        r["default_threshold"]
    )

    optimized_threshold = (
        r["optimized_threshold"]
    )

    uncalibrated = (
        r["uncalibrated"]
    )

    isotonic = (
        r["isotonic"]
    )

    final = r["final"]

    auc_difference = (
        epigenomic["roc_auc_mean"]
        - clinical["roc_auc_mean"]
    )

    results = f"""
# Results

## Predictive contribution of DNA methylation

Leakage-aware five-fold cross-validation demonstrated a substantial increase in discriminatory performance after incorporation of DNA methylation features. The clinical-genomic baseline model achieved a mean ROC AUC of {clinical["roc_auc_mean"]:.3f} ± {clinical["roc_auc_std"]:.3f}. In comparison, the clinical-genomic-epigenomic model achieved a mean ROC AUC of {epigenomic["roc_auc_mean"]:.3f} ± {epigenomic["roc_auc_std"]:.3f}. The absolute improvement in mean ROC AUC associated with inclusion of methylation features was therefore {auc_difference:.3f}.

The clinical-genomic-epigenomic model also achieved a mean balanced accuracy of {epigenomic["balanced_accuracy_mean"]:.3f} ± {epigenomic["balanced_accuracy_std"]:.3f}. These results indicated that the methylation layer contributed predictive information beyond that available from the clinical-genomic baseline.

## Effect of classification-threshold optimisation

At the default classification threshold of 0.50, the epigenomic model achieved a mean accuracy of {default_threshold["accuracy_mean"]:.3f}, but minority-class sensitivity remained low. Mean recall at the default threshold was {default_threshold["recall_mean"]:.3f} ± {default_threshold["recall_std"]:.3f}, while mean balanced accuracy was {default_threshold["balanced_accuracy_mean"]:.3f} ± {default_threshold["balanced_accuracy_std"]:.3f} and mean F1-score was {default_threshold["f1_score_mean"]:.3f} ± {default_threshold["f1_score_std"]:.3f}.

Threshold optimisation performed within the training folds identified a mean operating threshold of {optimized_threshold["threshold_mean"]:.3f} ± {optimized_threshold["threshold_std"]:.3f}. At this threshold, mean recall increased to {optimized_threshold["recall_mean"]:.3f} ± {optimized_threshold["recall_std"]:.3f}, mean balanced accuracy increased to {optimized_threshold["balanced_accuracy_mean"]:.3f} ± {optimized_threshold["balanced_accuracy_std"]:.3f}, and mean F1-score increased to {optimized_threshold["f1_score_mean"]:.3f} ± {optimized_threshold["f1_score_std"]:.3f}. Precision decreased from {default_threshold["precision_mean"]:.3f} at the default threshold to {optimized_threshold["precision_mean"]:.3f} after optimisation, demonstrating the expected trade-off between sensitivity and positive predictive value.

ROC AUC remained {optimized_threshold["roc_auc_mean"]:.3f}, and PR AUC remained {optimized_threshold["pr_auc_mean"]:.3f}, as these measures are not determined by a single classification threshold.

## Probability calibration

Probability calibration improved probabilistic model performance. The uncalibrated Random Forest produced a mean Brier score of {uncalibrated["brier_score_mean"]:.4f} ± {uncalibrated["brier_score_std"]:.4f} and a mean log loss of {uncalibrated["log_loss_mean"]:.4f} ± {uncalibrated["log_loss_std"]:.4f}.

Isotonic calibration produced the lowest mean Brier score, decreasing it to {isotonic["brier_score_mean"]:.4f} ± {isotonic["brier_score_std"]:.4f}. Mean log loss was also reduced to {isotonic["log_loss_mean"]:.4f} ± {isotonic["log_loss_std"]:.4f}. Following isotonic calibration, mean ROC AUC was {isotonic["roc_auc_mean"]:.3f} ± {isotonic["roc_auc_std"]:.3f} and mean PR AUC was {isotonic["pr_auc_mean"]:.3f} ± {isotonic["pr_auc_std"]:.3f}.

On the basis of the Brier score, isotonic regression was selected as the final probability-calibration method.

## Final development-stage model

The final model consisted of a Random Forest classifier combined with isotonic probability calibration. A calibrated operating threshold of {final["threshold"]:.3f} was selected using leakage-aware out-of-fold predictions.

The final out-of-fold ROC AUC was {final["roc_auc"]:.3f}, with a PR AUC of {final["pr_auc"]:.3f}. The final Brier score was {final["brier_score"]:.4f} and log loss was {final["log_loss"]:.4f}.

At the selected calibrated threshold, the model achieved an accuracy of {final["accuracy"]:.3f}, balanced accuracy of {final["balanced_accuracy"]:.3f}, precision of {final["precision"]:.3f}, recall of {final["recall"]:.3f}, and F1-score of {final["f1_score"]:.3f}.

Overall, the validation analyses showed that the initially observed near-perfect single-split performance was reduced under more stringent leakage-aware cross-validation, but the epigenomic model retained substantially stronger discrimination than the clinical-genomic baseline.
"""

    return results.strip()


# =====================================================
# DISCUSSION SECTION
# =====================================================

def create_discussion_section(
    r: dict,
) -> str:
    """Create manuscript-style Discussion section."""

    clinical = r["clinical"]
    epigenomic = r["epigenomic"]

    default_threshold = (
        r["default_threshold"]
    )

    optimized_threshold = (
        r["optimized_threshold"]
    )

    uncalibrated = (
        r["uncalibrated"]
    )

    isotonic = (
        r["isotonic"]
    )

    final = r["final"]

    auc_difference = (
        epigenomic["roc_auc_mean"]
        - clinical["roc_auc_mean"]
    )

    discussion = f"""
# Discussion

The integration of DNA methylation with clinical and genomic data resulted in a marked improvement in predictive discrimination. Under leakage-aware five-fold cross-validation, mean ROC AUC increased from {clinical["roc_auc_mean"]:.3f} for the clinical-genomic model to {epigenomic["roc_auc_mean"]:.3f} for the clinical-genomic-epigenomic model, corresponding to an absolute increase of {auc_difference:.3f}. This finding supports the interpretation that epigenomic variation contains information that is complementary to the baseline clinical and genomic feature set.

An important methodological observation was the difference between the original single-split model performance and the results obtained under stricter cross-validation. The earlier ROC AUC of approximately 1.0 was not reproduced when feature selection, preprocessing and model evaluation were constrained within cross-validation folds. The reduction to a mean ROC AUC of approximately {epigenomic["roc_auc_mean"]:.3f} illustrates the importance of leakage-aware validation in high-dimensional multi-omics analyses. Although the revised performance estimate was lower, it remained substantially above that of the clinical-genomic baseline.

Class imbalance also had a major effect on model interpretation. At the default probability threshold of 0.50, overall accuracy was high ({default_threshold["accuracy_mean"]:.3f}), but mean recall was only {default_threshold["recall_mean"]:.3f}. This demonstrated that conventional accuracy alone would have overstated model usefulness for detecting the minority outcome class. Optimising the operating threshold increased recall to {optimized_threshold["recall_mean"]:.3f} and balanced accuracy to {optimized_threshold["balanced_accuracy_mean"]:.3f}. The associated reduction in precision to {optimized_threshold["precision_mean"]:.3f} reflects the expected sensitivity-specificity trade-off and highlights the need to define the intended use case before selecting an operating threshold.

Probability calibration provided an additional improvement in model interpretability. Isotonic calibration reduced the mean Brier score from {uncalibrated["brier_score_mean"]:.4f} to {isotonic["brier_score_mean"]:.4f} and reduced mean log loss from {uncalibrated["log_loss_mean"]:.4f} to {isotonic["log_loss_mean"]:.4f}. The small reduction in ROC AUC following isotonic calibration is not inconsistent with successful calibration, because calibration primarily concerns agreement between predicted probabilities and observed event frequencies rather than ranking performance.

The final development-stage model achieved an out-of-fold ROC AUC of {final["roc_auc"]:.3f} and PR AUC of {final["pr_auc"]:.3f}, with a Brier score of {final["brier_score"]:.4f}. At the final calibrated threshold of {final["threshold"]:.3f}, balanced accuracy was {final["balanced_accuracy"]:.3f} and recall was {final["recall"]:.3f}. These values suggest that the model retained meaningful discrimination while providing a more balanced operating point than the original default-threshold classifier.

Several limitations remain. First, all estimates were obtained through internal validation within the same overall cohort. Internal cross-validation reduces optimism but does not establish transportability to an independent population. Second, the high dimensionality of methylation data relative to the number of patients increases the risk of overfitting and instability in individual selected features. Third, the operating threshold was selected to maximise balanced accuracy within the development data and should not be interpreted as a clinically validated cutoff. Fourth, calibration performance may change when disease prevalence or patient characteristics differ in external cohorts.

The results therefore support further investigation of methylation-enhanced multi-omics prediction, but they do not establish clinical utility. Independent external validation should be the next major analytical step. In addition, future work could examine stability of selected methylation loci, biological pathway enrichment, alternative feature-selection strategies and model architectures, and validation across independent breast cancer cohorts.
"""

    return discussion.strip()


# =====================================================
# CONCLUSIONS SECTION
# =====================================================

def create_conclusions_section(
    r: dict,
) -> str:
    """Create manuscript-style Conclusions section."""

    clinical = r["clinical"]
    epigenomic = r["epigenomic"]
    final = r["final"]

    conclusions = f"""
# Conclusions

Integration of DNA methylation into the clinical-genomic workflow substantially improved predictive discrimination in the available development cohort. Leakage-aware cross-validation increased mean ROC AUC from {clinical["roc_auc_mean"]:.3f} for the clinical-genomic baseline to {epigenomic["roc_auc_mean"]:.3f} after inclusion of epigenomic features.

Threshold optimisation improved detection of the minority outcome class, while isotonic calibration improved the quality of predicted probabilities. The final calibrated model achieved an internal out-of-fold ROC AUC of {final["roc_auc"]:.3f}, PR AUC of {final["pr_auc"]:.3f}, balanced accuracy of {final["balanced_accuracy"]:.3f}, and recall of {final["recall"]:.3f}.

These findings support the potential value of DNA methylation as an additional predictive layer in multi-omics modelling. However, the model remains a development-stage research model. Independent external validation is required before the reported probabilities, selected threshold, or model outputs can be considered generalisable or clinically applicable.
"""

    return conclusions.strip()


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Generate manuscript sections from Month 2 results."""

    print(
        "Loading Month 2 analytical results..."
    )

    results = load_all_results()

    print(
        "Generating Results section..."
    )

    results_text = (
        create_results_section(
            results
        )
    )

    RESULTS_FILE.write_text(
        results_text,
        encoding="utf-8",
    )

    print(
        "Generating Discussion section..."
    )

    discussion_text = (
        create_discussion_section(
            results
        )
    )

    DISCUSSION_FILE.write_text(
        discussion_text,
        encoding="utf-8",
    )

    print(
        "Generating Conclusions section..."
    )

    conclusions_text = (
        create_conclusions_section(
            results
        )
    )

    CONCLUSIONS_FILE.write_text(
        conclusions_text,
        encoding="utf-8",
    )

    combined = (
        results_text
        + "\n\n"
        + discussion_text
        + "\n\n"
        + conclusions_text
        + "\n"
    )

    COMBINED_FILE.write_text(
        combined,
        encoding="utf-8",
    )

    print(
        "\n========================================"
    )

    print(
        "MANUSCRIPT SECTIONS COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        "\nResults section saved to:\n"
        f"{RESULTS_FILE}"
    )

    print(
        "\nDiscussion section saved to:\n"
        f"{DISCUSSION_FILE}"
    )

    print(
        "\nConclusions section saved to:\n"
        f"{CONCLUSIONS_FILE}"
    )

    print(
        "\nCombined manuscript saved to:\n"
        f"{COMBINED_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()