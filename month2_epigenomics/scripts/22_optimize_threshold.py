from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pandas.api.types import is_numeric_dtype

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
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


PREDICTIONS_FILE = (
    TABLES_DIR
    / "threshold_validation_predictions.csv"
)

FOLD_RESULTS_FILE = (
    TABLES_DIR
    / "threshold_validation_folds.csv"
)

SUMMARY_FILE = (
    TABLES_DIR
    / "threshold_optimization_summary.csv"
)

THRESHOLD_FILE = (
    TABLES_DIR
    / "optimized_thresholds.csv"
)

ROC_FIGURE = (
    FIGURES_DIR
    / "threshold_validation_roc_curve.png"
)

PR_FIGURE = (
    FIGURES_DIR
    / "threshold_validation_precision_recall.png"
)

CONFUSION_DEFAULT_FIGURE = (
    FIGURES_DIR
    / "confusion_matrix_threshold_050.png"
)

CONFUSION_OPTIMAL_FIGURE = (
    FIGURES_DIR
    / "confusion_matrix_optimized_threshold.png"
)

REPORT_FILE = (
    REPORTS_DIR
    / "threshold_optimization_report.md"
)


# =====================================================
# SETTINGS
# =====================================================

PATIENT_COLUMN = "patient_id"
TARGET_COLUMN = "outcome"

OUTER_FOLDS = 5
INNER_FOLDS = 4

RANDOM_STATE = 42

TOP_METHYLATION_FEATURES = 500

N_ESTIMATORS = 300

MAX_CATEGORICAL_LEVELS = 30

MIN_VARIANCE = 1e-8

DEFAULT_THRESHOLD = 0.50


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
    """Add outcome from Month 1 if required."""

    if TARGET_COLUMN in epigenomic.columns:
        return epigenomic

    print(
        "Outcome not found in epigenomic table."
    )

    print(
        "Adding outcome from Month 1 table..."
    )

    baseline = standardize_patient_ids(
        baseline
    )

    epigenomic = standardize_patient_ids(
        epigenomic
    )

    if TARGET_COLUMN not in baseline.columns:
        raise KeyError(
            f"Month 1 table does not contain "
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
    """Clean patient IDs and binary outcome."""

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
    """Detect IDs and target-leaking variables."""

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

    # Survival leakage
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

    target_proxy_fragments = [
        "death",
        "deceased",
        "alive",
        "living",
        "mortality",
    ]

    if any(
        fragment in lower
        for fragment in target_proxy_fragments
    ):
        return True

    return False


def identify_safe_non_methylation_features(
    df: pd.DataFrame,
) -> list[str]:
    """Identify leakage-safe clinical/genomic features."""

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
    """Identify methylation features."""

    return [
        column
        for column in df.columns
        if str(column).startswith("meth_")
    ]


# =====================================================
# METHYLATION SELECTION
# =====================================================

def select_methylation_features(
    train_df: pd.DataFrame,
    methylation_columns: list[str],
) -> list[str]:
    """
    Select high-variance methylation features
    using training data only.
    """

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
# HIGH-CARDINALITY FILTER
# =====================================================

def identify_high_cardinality_columns(
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
# MODEL
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
            "No predictors available."
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
# THRESHOLD
# =====================================================

def find_best_threshold(
    y_true: pd.Series | np.ndarray,
    probabilities: np.ndarray,
) -> tuple[
    float,
    pd.DataFrame,
]:
    """
    Select threshold maximizing balanced accuracy.

    Threshold selection is done only with
    training-fold predictions.
    """

    thresholds = np.arange(
        0.05,
        0.951,
        0.01,
    )

    rows = []

    for threshold in thresholds:

        predictions = (
            probabilities
            >= threshold
        ).astype(int)

        balanced_accuracy = (
            balanced_accuracy_score(
                y_true,
                predictions,
            )
        )

        f1 = f1_score(
            y_true,
            predictions,
            zero_division=0,
        )

        recall = recall_score(
            y_true,
            predictions,
            zero_division=0,
        )

        precision = precision_score(
            y_true,
            predictions,
            zero_division=0,
        )

        rows.append(
            {
                "threshold":
                    threshold,
                "balanced_accuracy":
                    balanced_accuracy,
                "f1_score":
                    f1,
                "recall":
                    recall,
                "precision":
                    precision,
            }
        )

    table = pd.DataFrame(
        rows
    )

    best_row = (
        table
        .sort_values(
            [
                "balanced_accuracy",
                "f1_score",
            ],
            ascending=[
                False,
                False,
            ],
        )
        .iloc[0]
    )

    return (
        float(
            best_row[
                "threshold"
            ]
        ),
        table,
    )


# =====================================================
# METRICS
# =====================================================

def calculate_threshold_metrics(
    y_true: pd.Series | np.ndarray,
    probabilities: np.ndarray,
    threshold: float,
) -> dict[str, float]:
    """Calculate metrics for a threshold."""

    predictions = (
        probabilities
        >= threshold
    ).astype(int)

    return {
        "threshold":
            threshold,
        "accuracy":
            accuracy_score(
                y_true,
                predictions,
            ),
        "balanced_accuracy":
            balanced_accuracy_score(
                y_true,
                predictions,
            ),
        "precision":
            precision_score(
                y_true,
                predictions,
                zero_division=0,
            ),
        "recall":
            recall_score(
                y_true,
                predictions,
                zero_division=0,
            ),
        "f1_score":
            f1_score(
                y_true,
                predictions,
                zero_division=0,
            ),
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
    }


# =====================================================
# OUTER VALIDATION
# =====================================================

def run_nested_threshold_validation(
    dataset: pd.DataFrame,
    non_methylation_features: list[str],
    methylation_columns: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """
    Run outer CV with threshold tuning
    performed inside each training fold.
    """

    outer_cv = StratifiedKFold(
        n_splits=OUTER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    fold_results = []

    prediction_rows = []

    for fold_number, (
        train_index,
        test_index,
    ) in enumerate(
        outer_cv.split(
            dataset,
            dataset[
                TARGET_COLUMN
            ],
        ),
        start=1,
    ):

        print(
            "\n========================================"
        )

        print(
            f"OUTER FOLD {fold_number}/{OUTER_FOLDS}"
        )

        print(
            "========================================"
        )

        train_df = (
            dataset
            .iloc[
                train_index
            ]
            .reset_index(
                drop=True
            )
        )

        test_df = (
            dataset
            .iloc[
                test_index
            ]
            .reset_index(
                drop=True
            )
        )

        # ---------------------------------
        # FEATURE SELECTION
        # OUTER TRAIN ONLY
        # ---------------------------------

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

        # ---------------------------------
        # HIGH CARDINALITY
        # TRAIN ONLY
        # ---------------------------------

        high_cardinality = (
            identify_high_cardinality_columns(
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

        print(
            f"Selected methylation: "
            f"{len(selected_methylation):,}"
        )

        # ---------------------------------
        # INNER CV
        # THRESHOLD TUNING
        # ---------------------------------

        inner_cv = StratifiedKFold(
            n_splits=INNER_FOLDS,
            shuffle=True,
            random_state=(
                RANDOM_STATE
                + fold_number
            ),
        )

        threshold_pipeline = (
            build_pipeline(
                X_train
            )
        )

        print(
            "Generating inner-CV "
            "probabilities..."
        )

        inner_probabilities = (
            cross_val_predict(
                threshold_pipeline,
                X_train,
                y_train,
                cv=inner_cv,
                method="predict_proba",
                n_jobs=1,
            )[:, 1]
        )

        best_threshold, _ = (
            find_best_threshold(
                y_true=y_train,
                probabilities=(
                    inner_probabilities
                ),
            )
        )

        print(
            f"Optimal threshold from "
            f"training data: "
            f"{best_threshold:.2f}"
        )

        # ---------------------------------
        # TRAIN FINAL OUTER MODEL
        # ---------------------------------

        final_pipeline = (
            build_pipeline(
                X_train
            )
        )

        final_pipeline.fit(
            X_train,
            y_train,
        )

        test_probabilities = (
            final_pipeline
            .predict_proba(
                X_test
            )[:, 1]
        )

        # ---------------------------------
        # DEFAULT 0.50
        # ---------------------------------

        default_metrics = (
            calculate_threshold_metrics(
                y_true=y_test,
                probabilities=(
                    test_probabilities
                ),
                threshold=(
                    DEFAULT_THRESHOLD
                ),
            )
        )

        default_metrics.update(
            {
                "fold":
                    fold_number,
                "strategy":
                    "default_0.50",
            }
        )

        fold_results.append(
            default_metrics
        )

        # ---------------------------------
        # OPTIMIZED
        # ---------------------------------

        optimized_metrics = (
            calculate_threshold_metrics(
                y_true=y_test,
                probabilities=(
                    test_probabilities
                ),
                threshold=(
                    best_threshold
                ),
            )
        )

        optimized_metrics.update(
            {
                "fold":
                    fold_number,
                "strategy":
                    "optimized",
            }
        )

        fold_results.append(
            optimized_metrics
        )

        optimized_predictions = (
            test_probabilities
            >= best_threshold
        ).astype(int)

        default_predictions = (
            test_probabilities
            >= DEFAULT_THRESHOLD
        ).astype(int)

        for position in range(
            len(test_df)
        ):

            prediction_rows.append(
                {
                    "fold":
                        fold_number,
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
                    "probability":
                        float(
                            test_probabilities[
                                position
                            ]
                        ),
                    "optimized_threshold":
                        best_threshold,
                    "prediction_default":
                        int(
                            default_predictions[
                                position
                            ]
                        ),
                    "prediction_optimized":
                        int(
                            optimized_predictions[
                                position
                            ]
                        ),
                }
            )

        print(
            "\nDefault threshold 0.50:"
        )

        print(
            f"  Balanced accuracy: "
            f"{default_metrics['balanced_accuracy']:.3f}"
        )

        print(
            f"  Recall: "
            f"{default_metrics['recall']:.3f}"
        )

        print(
            f"  F1: "
            f"{default_metrics['f1_score']:.3f}"
        )

        print(
            "\nOptimized threshold:"
        )

        print(
            f"  Balanced accuracy: "
            f"{optimized_metrics['balanced_accuracy']:.3f}"
        )

        print(
            f"  Recall: "
            f"{optimized_metrics['recall']:.3f}"
        )

        print(
            f"  F1: "
            f"{optimized_metrics['f1_score']:.3f}"
        )

    return (
        pd.DataFrame(
            fold_results
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
    """Create mean and SD across outer folds."""

    metric_columns = [
        "accuracy",
        "balanced_accuracy",
        "precision",
        "recall",
        "f1_score",
        "roc_auc",
        "pr_auc",
        "threshold",
    ]

    rows = []

    for strategy in [
        "default_0.50",
        "optimized",
    ]:

        subset = fold_results[
            fold_results[
                "strategy"
            ] == strategy
        ]

        row = {
            "strategy":
                strategy,
            "folds":
                len(subset),
        }

        for metric in metric_columns:

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
# FIGURES
# =====================================================

def plot_roc_curve(
    predictions: pd.DataFrame,
) -> None:
    """Create pooled outer-fold ROC curve."""

    y_true = predictions[
        "true_outcome"
    ]

    probabilities = predictions[
        "probability"
    ]

    fpr, tpr, _ = roc_curve(
        y_true,
        probabilities,
    )

    roc_auc = roc_auc_score(
        y_true,
        probabilities,
    )

    plt.figure(
        figsize=(8, 6)
    )

    plt.plot(
        fpr,
        tpr,
        label=(
            f"Epigenomic model "
            f"(AUC = {roc_auc:.3f})"
        ),
    )

    plt.plot(
        [0, 1],
        [0, 1],
        linestyle="--",
        label="Random",
    )

    plt.xlabel(
        "False positive rate"
    )

    plt.ylabel(
        "True positive rate"
    )

    plt.title(
        "Leakage-aware outer-CV ROC curve"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        ROC_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_precision_recall_curve(
    predictions: pd.DataFrame,
) -> None:
    """Create pooled precision-recall curve."""

    y_true = predictions[
        "true_outcome"
    ]

    probabilities = predictions[
        "probability"
    ]

    precision, recall, _ = (
        precision_recall_curve(
            y_true,
            probabilities,
        )
    )

    pr_auc = (
        average_precision_score(
            y_true,
            probabilities,
        )
    )

    prevalence = (
        y_true.mean()
    )

    plt.figure(
        figsize=(8, 6)
    )

    plt.plot(
        recall,
        precision,
        label=(
            f"Epigenomic model "
            f"(AP = {pr_auc:.3f})"
        ),
    )

    plt.axhline(
        prevalence,
        linestyle="--",
        label=(
            f"Class prevalence "
            f"({prevalence:.3f})"
        ),
    )

    plt.xlabel(
        "Recall"
    )

    plt.ylabel(
        "Precision"
    )

    plt.title(
        "Leakage-aware precision-recall curve"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PR_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def save_confusion_matrices(
    predictions: pd.DataFrame,
) -> None:
    """Save default and optimized confusion matrices."""

    y_true = predictions[
        "true_outcome"
    ]

    default_predictions = predictions[
        "prediction_default"
    ]

    optimized_predictions = predictions[
        "prediction_optimized"
    ]

    default_matrix = confusion_matrix(
        y_true,
        default_predictions,
    )

    display = ConfusionMatrixDisplay(
        confusion_matrix=(
            default_matrix
        ),
    )

    display.plot()

    plt.title(
        "Confusion matrix - threshold 0.50"
    )

    plt.tight_layout()

    plt.savefig(
        CONFUSION_DEFAULT_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    optimized_matrix = confusion_matrix(
        y_true,
        optimized_predictions,
    )

    display = ConfusionMatrixDisplay(
        confusion_matrix=(
            optimized_matrix
        ),
    )

    display.plot()

    plt.title(
        "Confusion matrix - optimized threshold"
    )

    plt.tight_layout()

    plt.savefig(
        CONFUSION_OPTIMAL_FIGURE,
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
    """Safely convert dataframe to Markdown."""

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
    predictions: pd.DataFrame,
) -> str:
    """Create threshold-optimization report."""

    optimized_threshold = (
        fold_results.loc[
            fold_results[
                "strategy"
            ] == "optimized",
            "threshold",
        ]
        .mean()
    )

    pooled_roc_auc = (
        roc_auc_score(
            predictions[
                "true_outcome"
            ],
            predictions[
                "probability"
            ],
        )
    )

    pooled_pr_auc = (
        average_precision_score(
            predictions[
                "true_outcome"
            ],
            predictions[
                "probability"
            ],
        )
    )

    prevalence = (
        predictions[
            "true_outcome"
        ]
        .mean()
    )

    report_lines = [
        "# Classification Threshold Optimisation Report",
        "",
        "## Objective",
        "",
        (
            "The leakage-aware validation in Step 21 showed "
            "that the epigenomic model achieved strong ROC AUC "
            "but relatively low recall when the default "
            "classification threshold of 0.50 was used."
        ),
        "",
        (
            "This analysis evaluates whether a training-derived "
            "decision threshold improves classification of the "
            "minority outcome class."
        ),
        "",
        "## Validation design",
        "",
        (
            f"- Patients: "
            f"{len(dataset):,}"
        ),
        (
            f"- Outer validation folds: "
            f"{OUTER_FOLDS}"
        ),
        (
            f"- Inner threshold-tuning folds: "
            f"{INNER_FOLDS}"
        ),
        (
            "- Methylation feature selection was performed "
            "using outer-training data only."
        ),
        (
            "- Threshold optimisation was performed using "
            "predictions generated only from training data."
        ),
        (
            "- The optimized threshold was then evaluated "
            "on the independent outer test fold."
        ),
        "",
        "## Outcome distribution",
        "",
        (
            f"- Positive-class prevalence: "
            f"{prevalence:.3f}"
        ),
        "",
        "## Cross-validation results",
        "",
        dataframe_to_markdown(
            summary.round(4)
        ),
        "",
        "## Fold-level results",
        "",
        dataframe_to_markdown(
            fold_results[
                [
                    "fold",
                    "strategy",
                    "threshold",
                    "accuracy",
                    "balanced_accuracy",
                    "precision",
                    "recall",
                    "f1_score",
                    "roc_auc",
                    "pr_auc",
                ]
            ].round(4)
        ),
        "",
        "## Probability discrimination",
        "",
        (
            f"- Pooled outer-fold ROC AUC: "
            f"{pooled_roc_auc:.3f}"
        ),
        (
            f"- Pooled outer-fold PR AUC: "
            f"{pooled_pr_auc:.3f}"
        ),
        "",
        "## Optimized threshold",
        "",
        (
            "- Mean training-derived optimized threshold: "
            f"{optimized_threshold:.3f}"
        ),
        "",
        "## Interpretation",
        "",
        (
            "ROC AUC describes ranking performance and does not "
            "depend on a single classification threshold."
        ),
        "",
        (
            "Because the outcome classes are imbalanced, accuracy "
            "alone can be misleading. Balanced accuracy, recall, "
            "F1-score and precision-recall performance should "
            "therefore be considered alongside ROC AUC."
        ),
        "",
        (
            "The optimized threshold should only be preferred "
            "if it improves clinically relevant metrics such as "
            "recall and balanced accuracy on the independent "
            "outer folds."
        ),
        "",
        (
            "The threshold estimated here should be considered "
            "a validation-stage threshold rather than a final "
            "clinical decision threshold."
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
    """Optimize classification threshold safely."""

    print(
        "Loading Month 1 outcome table..."
    )

    baseline = load_dataset(
        BASELINE_FILE
    )

    print(
        "Loading integrated epigenomic data..."
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
        identify_safe_non_methylation_features(
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
            "No methylation features were found."
        )

    # =================================================
    # NESTED THRESHOLD VALIDATION
    # =================================================

    (
        fold_results,
        predictions,
    ) = run_nested_threshold_validation(
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
    # SAVE TABLES
    # =================================================

    fold_results.to_csv(
        FOLD_RESULTS_FILE,
        index=False,
    )

    predictions.to_csv(
        PREDICTIONS_FILE,
        index=False,
    )

    summary.to_csv(
        SUMMARY_FILE,
        index=False,
    )

    optimized_thresholds = (
        fold_results[
            fold_results[
                "strategy"
            ] == "optimized"
        ][
            [
                "fold",
                "threshold",
                "balanced_accuracy",
                "precision",
                "recall",
                "f1_score",
            ]
        ]
        .copy()
    )

    optimized_thresholds.to_csv(
        THRESHOLD_FILE,
        index=False,
    )

    # =================================================
    # FIGURES
    # =================================================

    plot_roc_curve(
        predictions
    )

    plot_precision_recall_curve(
        predictions
    )

    save_confusion_matrices(
        predictions
    )

    # =================================================
    # REPORT
    # =================================================

    report = create_report(
        dataset=dataset,
        fold_results=fold_results,
        summary=summary,
        predictions=predictions,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # =================================================
    # TERMINAL OUTPUT
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "THRESHOLD OPTIMIZATION SUMMARY"
    )

    print(
        "========================================"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    mean_threshold = (
        optimized_thresholds[
            "threshold"
        ]
        .mean()
    )

    print(
        "\nMean optimized threshold: "
        f"{mean_threshold:.3f}"
    )

    print(
        "\nFinished successfully."
    )

    print(
        "\nFold results saved to:\n"
        f"{FOLD_RESULTS_FILE}"
    )

    print(
        "\nPredictions saved to:\n"
        f"{PREDICTIONS_FILE}"
    )

    print(
        "\nSummary saved to:\n"
        f"{SUMMARY_FILE}"
    )

    print(
        "\nOptimized thresholds saved to:\n"
        f"{THRESHOLD_FILE}"
    )

    print(
        "\nROC curve saved to:\n"
        f"{ROC_FIGURE}"
    )

    print(
        "\nPrecision-recall curve saved to:\n"
        f"{PR_FIGURE}"
    )

    print(
        "\nDefault confusion matrix saved to:\n"
        f"{CONFUSION_DEFAULT_FIGURE}"
    )

    print(
        "\nOptimized confusion matrix saved to:\n"
        f"{CONFUSION_OPTIMAL_FIGURE}"
    )

    print(
        "\nReport saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()