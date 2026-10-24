from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pandas.api.types import is_numeric_dtype

from sklearn.calibration import (
    CalibratedClassifierCV,
    calibration_curve,
)
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)
from sklearn.model_selection import (
    StratifiedKFold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
)

TABLES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "tables"
)

FIGURES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "figures"
)

REPORTS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
)

for directory in [
    TABLES_DIR,
    FIGURES_DIR,
    REPORTS_DIR,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


BASELINE_FILE = (
    PROCESSED_DIR
    / "ml_feature_table.csv"
)

EPIGENOMIC_FILE = (
    PROCESSED_DIR
    / "clinical_genomic_epigenomic_merged.csv"
)

FOLD_RESULTS_FILE = (
    TABLES_DIR
    / "calibration_validation_folds.csv"
)

SUMMARY_FILE = (
    TABLES_DIR
    / "calibration_validation_summary.csv"
)

PREDICTIONS_FILE = (
    TABLES_DIR
    / "calibration_predictions.csv"
)

CALIBRATION_CURVE_FILE = (
    FIGURES_DIR
    / "model_calibration_curve.png"
)

PROBABILITY_HISTOGRAM_FILE = (
    FIGURES_DIR
    / "calibrated_probability_distribution.png"
)

REPORT_FILE = (
    REPORTS_DIR
    / "model_calibration_report.md"
)


# =====================================================
# SETTINGS
# =====================================================

PATIENT_COLUMN = "patient_id"
TARGET_COLUMN = "outcome"

OUTER_FOLDS = 5
CALIBRATION_FOLDS = 4

RANDOM_STATE = 42

TOP_METHYLATION_FEATURES = 500

N_ESTIMATORS = 300

MAX_CATEGORICAL_LEVELS = 30

MIN_VARIANCE = 1e-8

N_CALIBRATION_BINS = 10


# =====================================================
# DATA LOADING
# =====================================================

def load_dataset(
    file_path: Path,
) -> pd.DataFrame:
    """Load CSV dataset."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def standardize_patient_ids(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Standardize patient identifiers."""

    if PATIENT_COLUMN not in df.columns:
        raise KeyError(
            f"Missing column: {PATIENT_COLUMN}"
        )

    df = df.copy()

    df[PATIENT_COLUMN] = (
        df[PATIENT_COLUMN]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return df


# =====================================================
# OUTCOME
# =====================================================

def add_outcome(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
) -> pd.DataFrame:
    """Add outcome from Month 1 if missing."""

    if TARGET_COLUMN in epigenomic.columns:
        return epigenomic

    print(
        "Outcome not found in epigenomic dataset."
    )

    print(
        "Adding outcome from Month 1 dataset..."
    )

    baseline = standardize_patient_ids(
        baseline
    )

    epigenomic = standardize_patient_ids(
        epigenomic
    )

    if TARGET_COLUMN not in baseline.columns:
        raise KeyError(
            f"Month 1 dataset does not contain "
            f"'{TARGET_COLUMN}'."
        )

    outcome_table = baseline[
        [
            PATIENT_COLUMN,
            TARGET_COLUMN,
        ]
    ].copy()

    outcome_table[TARGET_COLUMN] = (
        pd.to_numeric(
            outcome_table[TARGET_COLUMN],
            errors="coerce",
        )
    )

    outcome_table = (
        outcome_table
        .dropna(
            subset=[
                PATIENT_COLUMN,
                TARGET_COLUMN,
            ]
        )
        .drop_duplicates(
            subset=[PATIENT_COLUMN],
            keep="first",
        )
    )

    epigenomic = epigenomic.merge(
        outcome_table,
        on=PATIENT_COLUMN,
        how="left",
        validate="one_to_one",
    )

    print(
        "Missing outcomes after matching: "
        f"{epigenomic[TARGET_COLUMN].isna().sum():,}"
    )

    return epigenomic


def clean_dataset(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Clean binary outcome and patient IDs."""

    df = standardize_patient_ids(
        df
    )

    df[TARGET_COLUMN] = (
        pd.to_numeric(
            df[TARGET_COLUMN],
            errors="coerce",
        )
    )

    df = df.dropna(
        subset=[
            PATIENT_COLUMN,
            TARGET_COLUMN,
        ]
    ).copy()

    df[TARGET_COLUMN] = (
        df[TARGET_COLUMN]
        .astype(int)
    )

    df = df[
        df[TARGET_COLUMN].isin(
            [0, 1]
        )
    ].copy()

    df = (
        df.drop_duplicates(
            subset=[PATIENT_COLUMN],
            keep="first",
        )
        .reset_index(drop=True)
    )

    return df


# =====================================================
# LEAKAGE CONTROL
# =====================================================

def is_forbidden_feature(
    column: str,
) -> bool:
    """Identify identifiers and potential target leakage."""

    name = str(column)
    lower = name.lower()

    exact_excluded = {
        PATIENT_COLUMN,
        TARGET_COLUMN,
        "sample_id",
        "case_id",
        "study_id",
        "project_id",
        "submitter_id",
        "patient_barcode",
        "sample_barcode",
        "vital_status",
        "survival_status",
        "os_status",
        "os_status_patient",
    }

    if name in exact_excluded:
        return True

    identifier_fragments = [
        "patient_id",
        "sample_id",
        "case_id",
        "barcode",
        "submitter_id",
        "project_id",
        "study_id",
        "tissue_source_site",
    ]

    if any(
        fragment in lower
        for fragment in identifier_fragments
    ):
        return True

    if lower.endswith("_months"):
        return True

    if "survival" in lower:
        return True

    if lower.startswith("os_"):
        return True

    if lower.startswith("dfs_"):
        return True

    if lower.startswith("pfs_"):
        return True

    if lower.startswith("dss_"):
        return True

    target_proxies = [
        "death",
        "deceased",
        "alive",
        "living",
        "mortality",
    ]

    if any(
        fragment in lower
        for fragment in target_proxies
    ):
        return True

    return False


def identify_non_methylation_features(
    df: pd.DataFrame,
) -> list[str]:
    """Identify safe non-methylation predictors."""

    return [
        column
        for column in df.columns
        if (
            not str(column).startswith("meth_")
            and not is_forbidden_feature(column)
        )
    ]


def identify_methylation_columns(
    df: pd.DataFrame,
) -> list[str]:
    """Identify methylation predictors."""

    return [
        column
        for column in df.columns
        if str(column).startswith("meth_")
    ]


# =====================================================
# METHYLATION FEATURE SELECTION
# =====================================================

def select_methylation_features(
    train_df: pd.DataFrame,
    methylation_columns: list[str],
) -> list[str]:
    """Select methylation features using training fold only."""

    methylation = train_df[
        methylation_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    variance = methylation.var(
        axis=0,
        skipna=True,
    )

    variance = variance[
        variance.notna()
        & (variance > MIN_VARIANCE)
    ]

    return (
        variance
        .sort_values(
            ascending=False
        )
        .head(
            TOP_METHYLATION_FEATURES
        )
        .index
        .tolist()
    )


# =====================================================
# CARDINALITY FILTER
# =====================================================

def identify_high_cardinality(
    train_df: pd.DataFrame,
) -> list[str]:
    """Identify categorical variables with too many levels."""

    categorical = [
        column
        for column in train_df.columns
        if not is_numeric_dtype(
            train_df[column]
        )
    ]

    return [
        column
        for column in categorical
        if train_df[column].nunique(
            dropna=True
        ) > MAX_CATEGORICAL_LEVELS
    ]


# =====================================================
# MODEL PIPELINE
# =====================================================

def build_pipeline(
    X_train: pd.DataFrame,
) -> Pipeline:
    """Build Random Forest preprocessing pipeline."""

    numeric_features = [
        column
        for column in X_train.columns
        if is_numeric_dtype(
            X_train[column]
        )
    ]

    categorical_features = [
        column
        for column in X_train.columns
        if column not in numeric_features
    ]

    transformers = []

    if numeric_features:

        numeric_pipeline = Pipeline(
            steps=[
                (
                    "imputer",
                    SimpleImputer(
                        strategy="median",
                        keep_empty_features=True,
                    ),
                ),
            ]
        )

        transformers.append(
            (
                "numeric",
                numeric_pipeline,
                numeric_features,
            )
        )

    if categorical_features:

        categorical_pipeline = Pipeline(
            steps=[
                (
                    "imputer",
                    SimpleImputer(
                        strategy="most_frequent",
                        keep_empty_features=True,
                    ),
                ),
                (
                    "encoder",
                    OneHotEncoder(
                        handle_unknown="ignore",
                    ),
                ),
            ]
        )

        transformers.append(
            (
                "categorical",
                categorical_pipeline,
                categorical_features,
            )
        )

    if not transformers:
        raise ValueError(
            "No model predictors remain."
        )

    preprocessor = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
    )

    classifier = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )

    return Pipeline(
        steps=[
            (
                "preprocessor",
                preprocessor,
            ),
            (
                "classifier",
                classifier,
            ),
        ]
    )


# =====================================================
# METRICS
# =====================================================

def calculate_probability_metrics(
    y_true,
    probabilities,
) -> dict[str, float]:
    """Calculate probability-based model metrics."""

    clipped = np.clip(
        probabilities,
        1e-7,
        1 - 1e-7,
    )

    return {
        "roc_auc":
            roc_auc_score(
                y_true,
                probabilities,
            ),
        "pr_auc":
            average_precision_score(
                y_true,
                probabilities,
            ),
        "brier_score":
            brier_score_loss(
                y_true,
                probabilities,
            ),
        "log_loss":
            log_loss(
                y_true,
                clipped,
                labels=[0, 1],
            ),
    }


# =====================================================
# CALIBRATION VALIDATION
# =====================================================

def run_calibration_validation(
    dataset: pd.DataFrame,
    non_methylation_features: list[str],
    methylation_columns: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Run leakage-aware outer CV calibration comparison."""

    outer_cv = StratifiedKFold(
        n_splits=OUTER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    fold_rows = []
    prediction_rows = []

    for fold, (
        train_index,
        test_index,
    ) in enumerate(
        outer_cv.split(
            dataset,
            dataset[TARGET_COLUMN],
        ),
        start=1,
    ):

        print(
            "\n========================================"
        )

        print(
            f"CALIBRATION FOLD "
            f"{fold}/{OUTER_FOLDS}"
        )

        print(
            "========================================"
        )

        train_df = (
            dataset
            .iloc[train_index]
            .reset_index(drop=True)
        )

        test_df = (
            dataset
            .iloc[test_index]
            .reset_index(drop=True)
        )

        selected_methylation = (
            select_methylation_features(
                train_df,
                methylation_columns,
            )
        )

        feature_columns = (
            non_methylation_features
            + selected_methylation
        )

        feature_columns = list(
            dict.fromkeys(
                feature_columns
            )
        )

        X_train = train_df[
            feature_columns
        ].copy()

        X_test = test_df[
            feature_columns
        ].copy()

        y_train = train_df[
            TARGET_COLUMN
        ].copy()

        y_test = test_df[
            TARGET_COLUMN
        ].copy()

        high_cardinality = (
            identify_high_cardinality(
                X_train
            )
        )

        if high_cardinality:

            X_train = X_train.drop(
                columns=high_cardinality
            )

            X_test = X_test.drop(
                columns=high_cardinality
            )

        print(
            f"Training patients: "
            f"{len(X_train):,}"
        )

        print(
            f"Testing patients: "
            f"{len(X_test):,}"
        )

        print(
            f"Predictors: "
            f"{X_train.shape[1]:,}"
        )

        # ---------------------------------
        # RAW RANDOM FOREST
        # ---------------------------------

        raw_model = build_pipeline(
            X_train
        )

        raw_model.fit(
            X_train,
            y_train,
        )

        raw_probabilities = (
            raw_model
            .predict_proba(
                X_test
            )[:, 1]
        )

        raw_metrics = (
            calculate_probability_metrics(
                y_test,
                raw_probabilities,
            )
        )

        raw_metrics.update(
            {
                "fold":
                    fold,
                "method":
                    "uncalibrated",
            }
        )

        fold_rows.append(
            raw_metrics
        )

        # ---------------------------------
        # SIGMOID CALIBRATION
        # ---------------------------------

        sigmoid_base = (
            build_pipeline(
                X_train
            )
        )

        sigmoid_model = (
            CalibratedClassifierCV(
                estimator=sigmoid_base,
                method="sigmoid",
                cv=CALIBRATION_FOLDS,
                n_jobs=1,
            )
        )

        sigmoid_model.fit(
            X_train,
            y_train,
        )

        sigmoid_probabilities = (
            sigmoid_model
            .predict_proba(
                X_test
            )[:, 1]
        )

        sigmoid_metrics = (
            calculate_probability_metrics(
                y_test,
                sigmoid_probabilities,
            )
        )

        sigmoid_metrics.update(
            {
                "fold":
                    fold,
                "method":
                    "sigmoid",
            }
        )

        fold_rows.append(
            sigmoid_metrics
        )

        # ---------------------------------
        # ISOTONIC CALIBRATION
        # ---------------------------------

        isotonic_base = (
            build_pipeline(
                X_train
            )
        )

        isotonic_model = (
            CalibratedClassifierCV(
                estimator=isotonic_base,
                method="isotonic",
                cv=CALIBRATION_FOLDS,
                n_jobs=1,
            )
        )

        isotonic_model.fit(
            X_train,
            y_train,
        )

        isotonic_probabilities = (
            isotonic_model
            .predict_proba(
                X_test
            )[:, 1]
        )

        isotonic_metrics = (
            calculate_probability_metrics(
                y_test,
                isotonic_probabilities,
            )
        )

        isotonic_metrics.update(
            {
                "fold":
                    fold,
                "method":
                    "isotonic",
            }
        )

        fold_rows.append(
            isotonic_metrics
        )

        # ---------------------------------
        # SAVE PATIENT PREDICTIONS
        # ---------------------------------

        for position in range(
            len(test_df)
        ):

            prediction_rows.append(
                {
                    "fold":
                        fold,
                    PATIENT_COLUMN:
                        test_df.iloc[
                            position
                        ][PATIENT_COLUMN],
                    "true_outcome":
                        int(
                            y_test.iloc[
                                position
                            ]
                        ),
                    "uncalibrated_probability":
                        float(
                            raw_probabilities[
                                position
                            ]
                        ),
                    "sigmoid_probability":
                        float(
                            sigmoid_probabilities[
                                position
                            ]
                        ),
                    "isotonic_probability":
                        float(
                            isotonic_probabilities[
                                position
                            ]
                        ),
                }
            )

        print(
            "\nBrier scores:"
        )

        print(
            f"  Uncalibrated: "
            f"{raw_metrics['brier_score']:.4f}"
        )

        print(
            f"  Sigmoid:     "
            f"{sigmoid_metrics['brier_score']:.4f}"
        )

        print(
            f"  Isotonic:    "
            f"{isotonic_metrics['brier_score']:.4f}"
        )

    return (
        pd.DataFrame(
            fold_rows
        ),
        pd.DataFrame(
            prediction_rows
        ),
    )


# =====================================================
# SUMMARY
# =====================================================

def create_summary(
    fold_results: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize calibration performance."""

    metrics = [
        "roc_auc",
        "pr_auc",
        "brier_score",
        "log_loss",
    ]

    rows = []

    for method in [
        "uncalibrated",
        "sigmoid",
        "isotonic",
    ]:

        subset = fold_results[
            fold_results[
                "method"
            ] == method
        ]

        row = {
            "method":
                method,
            "folds":
                len(subset),
        }

        for metric in metrics:

            row[
                f"{metric}_mean"
            ] = (
                subset[
                    metric
                ].mean()
            )

            row[
                f"{metric}_std"
            ] = (
                subset[
                    metric
                ].std()
            )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )


# =====================================================
# CALIBRATION CURVE
# =====================================================

def plot_calibration_curves(
    predictions: pd.DataFrame,
) -> None:
    """Plot pooled outer-fold calibration curves."""

    y_true = predictions[
        "true_outcome"
    ]

    probability_columns = {
        "Uncalibrated":
            "uncalibrated_probability",
        "Sigmoid":
            "sigmoid_probability",
        "Isotonic":
            "isotonic_probability",
    }

    plt.figure(
        figsize=(8, 7)
    )

    plt.plot(
        [0, 1],
        [0, 1],
        linestyle="--",
        label="Perfect calibration",
    )

    for label, column in (
        probability_columns.items()
    ):

        probability_true, (
            probability_predicted
        ) = calibration_curve(
            y_true,
            predictions[column],
            n_bins=N_CALIBRATION_BINS,
            strategy="quantile",
        )

        plt.plot(
            probability_predicted,
            probability_true,
            marker="o",
            label=label,
        )

    plt.xlabel(
        "Mean predicted probability"
    )

    plt.ylabel(
        "Observed outcome frequency"
    )

    plt.title(
        "Probability calibration curve"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        CALIBRATION_CURVE_FILE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_probability_histogram(
    predictions: pd.DataFrame,
) -> None:
    """Plot calibrated probability distributions."""

    plt.figure(
        figsize=(9, 6)
    )

    plt.hist(
        predictions[
            "uncalibrated_probability"
        ],
        bins=30,
        alpha=0.5,
        label="Uncalibrated",
    )

    plt.hist(
        predictions[
            "sigmoid_probability"
        ],
        bins=30,
        alpha=0.5,
        label="Sigmoid",
    )

    plt.hist(
        predictions[
            "isotonic_probability"
        ],
        bins=30,
        alpha=0.5,
        label="Isotonic",
    )

    plt.xlabel(
        "Predicted probability"
    )

    plt.ylabel(
        "Patients"
    )

    plt.title(
        "Predicted probability distributions"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PROBABILITY_HISTOGRAM_FILE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# =====================================================
# REPORT
# =====================================================

def dataframe_to_markdown(
    df: pd.DataFrame,
) -> str:
    """Convert dataframe to Markdown safely."""

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


def create_report(
    dataset: pd.DataFrame,
    fold_results: pd.DataFrame,
    summary: pd.DataFrame,
) -> str:
    """Create calibration report."""

    best_brier_row = (
        summary
        .sort_values(
            "brier_score_mean",
            ascending=True,
        )
        .iloc[0]
    )

    best_logloss_row = (
        summary
        .sort_values(
            "log_loss_mean",
            ascending=True,
        )
        .iloc[0]
    )

    report_lines = [
        "# Model Probability Calibration Report",
        "",
        "## Objective",
        "",
        (
            "This analysis evaluates whether probability "
            "calibration improves the reliability of predicted "
            "outcome probabilities produced by the "
            "clinical-genomic-epigenomic Random Forest model."
        ),
        "",
        "## Cohort",
        "",
        (
            f"- Patients: "
            f"{len(dataset):,}"
        ),
        (
            f"- Positive outcomes: "
            f"{int(dataset[TARGET_COLUMN].sum()):,}"
        ),
        (
            f"- Negative outcomes: "
            f"{int((dataset[TARGET_COLUMN] == 0).sum()):,}"
        ),
        "",
        "## Validation design",
        "",
        (
            f"- Outer validation folds: "
            f"{OUTER_FOLDS}"
        ),
        (
            f"- Calibration folds: "
            f"{CALIBRATION_FOLDS}"
        ),
        (
            "- Methylation feature selection was "
            "performed using outer-training data only."
        ),
        (
            "- Calibration was fitted without access "
            "to outer-test outcomes."
        ),
        "",
        "## Calibration methods",
        "",
        "- Uncalibrated Random Forest",
        "- Sigmoid / Platt-style calibration",
        "- Isotonic regression calibration",
        "",
        "## Cross-validation summary",
        "",
        dataframe_to_markdown(
            summary.round(4)
        ),
        "",
        "## Fold-level results",
        "",
        dataframe_to_markdown(
            fold_results.round(4)
        ),
        "",
        "## Best calibration results",
        "",
        (
            "- Lowest mean Brier score: "
            f"`{best_brier_row['method']}` "
            f"({best_brier_row['brier_score_mean']:.4f})"
        ),
        (
            "- Lowest mean log loss: "
            f"`{best_logloss_row['method']}` "
            f"({best_logloss_row['log_loss_mean']:.4f})"
        ),
        "",
        "## Interpretation",
        "",
        (
            "Lower Brier score and lower log loss indicate "
            "better probabilistic predictions."
        ),
        "",
        (
            "ROC AUC and PR AUC primarily describe discrimination, "
            "whereas calibration assesses whether predicted "
            "probabilities correspond to observed outcome frequencies."
        ),
        "",
        (
            "A calibrated probability should not automatically "
            "be interpreted as an externally validated clinical "
            "risk estimate. Independent validation would still "
            "be required before clinical use."
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
    """Validate probability calibration."""

    print(
        "Loading Month 1 outcome dataset..."
    )

    baseline = load_dataset(
        BASELINE_FILE
    )

    print(
        "Loading epigenomic dataset..."
    )

    dataset = load_dataset(
        EPIGENOMIC_FILE
    )

    baseline = clean_dataset(
        baseline
    )

    dataset = add_outcome(
        baseline=baseline,
        epigenomic=dataset,
    )

    dataset = clean_dataset(
        dataset
    )

    print(
        f"\nPatients available: "
        f"{len(dataset):,}"
    )

    print(
        "\nOutcome distribution:"
    )

    print(
        dataset[
            TARGET_COLUMN
        ]
        .value_counts()
        .sort_index()
    )

    non_methylation_features = (
        identify_non_methylation_features(
            dataset
        )
    )

    methylation_columns = (
        identify_methylation_columns(
            dataset
        )
    )

    print(
        f"\nSafe clinical/genomic features: "
        f"{len(non_methylation_features):,}"
    )

    print(
        f"Available methylation features: "
        f"{len(methylation_columns):,}"
    )

    if not methylation_columns:
        raise ValueError(
            "No methylation features found."
        )

    # =================================================
    # VALIDATION
    # =================================================

    (
        fold_results,
        predictions,
    ) = run_calibration_validation(
        dataset=dataset,
        non_methylation_features=(
            non_methylation_features
        ),
        methylation_columns=(
            methylation_columns
        ),
    )

    # =================================================
    # SUMMARY
    # =================================================

    summary = create_summary(
        fold_results
    )

    # =================================================
    # SAVE
    # =================================================

    fold_results.to_csv(
        FOLD_RESULTS_FILE,
        index=False,
    )

    summary.to_csv(
        SUMMARY_FILE,
        index=False,
    )

    predictions.to_csv(
        PREDICTIONS_FILE,
        index=False,
    )

    # =================================================
    # FIGURES
    # =================================================

    plot_calibration_curves(
        predictions
    )

    plot_probability_histogram(
        predictions
    )

    # =================================================
    # REPORT
    # =================================================

    report = create_report(
        dataset=dataset,
        fold_results=fold_results,
        summary=summary,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # =================================================
    # TERMINAL
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "CALIBRATION SUMMARY"
    )

    print(
        "========================================"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    best_brier = (
        summary
        .sort_values(
            "brier_score_mean"
        )
        .iloc[0]
    )

    print(
        "\nBest method by Brier score: "
        f"{best_brier['method']}"
    )

    print(
        "\nFinished successfully."
    )

    print(
        "\nFold results saved to:\n"
        f"{FOLD_RESULTS_FILE}"
    )

    print(
        "\nSummary saved to:\n"
        f"{SUMMARY_FILE}"
    )

    print(
        "\nPredictions saved to:\n"
        f"{PREDICTIONS_FILE}"
    )

    print(
        "\nCalibration curve saved to:\n"
        f"{CALIBRATION_CURVE_FILE}"
    )

    print(
        "\nProbability histogram saved to:\n"
        f"{PROBABILITY_HISTOGRAM_FILE}"
    )

    print(
        "\nReport saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()