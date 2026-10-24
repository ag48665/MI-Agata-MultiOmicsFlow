from pathlib import Path
import json

import joblib
import numpy as np
import pandas as pd

from pandas.api.types import is_numeric_dtype

from sklearn.calibration import CalibratedClassifierCV
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    log_loss,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
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

OUTPUTS_DIR = (
    PROJECT_ROOT
    / "outputs"
)

MODELS_DIR = (
    OUTPUTS_DIR
    / "models"
)

TABLES_DIR = (
    OUTPUTS_DIR
    / "tables"
)

REPORTS_DIR = (
    OUTPUTS_DIR
    / "reports"
)

for directory in [
    MODELS_DIR,
    TABLES_DIR,
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


FINAL_MODEL_FILE = (
    MODELS_DIR
    / "final_epigenomic_calibrated_model.pkl"
)

FINAL_CONFIG_FILE = (
    MODELS_DIR
    / "final_epigenomic_model_config.json"
)

FINAL_METRICS_FILE = (
    TABLES_DIR
    / "final_epigenomic_model_metrics.csv"
)

FINAL_FEATURES_FILE = (
    TABLES_DIR
    / "final_epigenomic_model_features.csv"
)

FINAL_THRESHOLD_FILE = (
    TABLES_DIR
    / "final_calibrated_threshold.csv"
)

OOF_PREDICTIONS_FILE = (
    TABLES_DIR
    / "final_calibrated_oof_predictions.csv"
)

THRESHOLD_SEARCH_FILE = (
    TABLES_DIR
    / "final_calibrated_threshold_search.csv"
)

REPORT_FILE = (
    REPORTS_DIR
    / "final_epigenomic_model_report.md"
)


# =====================================================
# SETTINGS
# =====================================================

PATIENT_COLUMN = "patient_id"
TARGET_COLUMN = "outcome"

RANDOM_STATE = 42

OUTER_FOLDS = 5
CALIBRATION_FOLDS = 4

TOP_METHYLATION_FEATURES = 500

N_ESTIMATORS = 300

MAX_CATEGORICAL_LEVELS = 30

MIN_VARIANCE = 1e-8

CALIBRATION_METHOD = "isotonic"


# =====================================================
# LOAD DATA
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
    """Standardize TCGA patient IDs."""

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
    """Add outcome from Month 1 if not present."""

    if TARGET_COLUMN in epigenomic.columns:
        print(
            "Outcome column already exists."
        )
        return epigenomic

    print(
        "Outcome column not found."
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
            f"'{TARGET_COLUMN}' not found "
            "in ml_feature_table.csv."
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

    missing = int(
        epigenomic[
            TARGET_COLUMN
        ]
        .isna()
        .sum()
    )

    print(
        f"Missing outcomes after merge: "
        f"{missing:,}"
    )

    return epigenomic


def clean_dataset(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Clean patient IDs and binary outcome."""

    df = standardize_patient_ids(
        df
    )

    if TARGET_COLUMN not in df.columns:
        raise KeyError(
            f"Missing target column: "
            f"{TARGET_COLUMN}"
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
    """Identify IDs and potential outcome leakage."""

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

    # Survival / follow-up leakage
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

    # Direct target proxies
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


def identify_non_methylation_features(
    df: pd.DataFrame,
) -> list[str]:
    """Identify safe clinical/genomic predictors."""

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
    """Select highest-variance methylation features."""

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

    selected = (
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

    return selected


# =====================================================
# HIGH-CARDINALITY CONTROL
# =====================================================

def identify_high_cardinality_columns(
    train_df: pd.DataFrame,
) -> list[str]:
    """Identify categorical predictors with too many levels."""

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

def build_base_pipeline(
    X_train: pd.DataFrame,
) -> Pipeline:
    """Create Random Forest preprocessing pipeline."""

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
            "No predictors remain after QC."
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


def build_calibrated_model(
    X_train: pd.DataFrame,
) -> CalibratedClassifierCV:
    """Build isotonic calibrated Random Forest."""

    base_model = build_base_pipeline(
        X_train
    )

    calibration_cv = StratifiedKFold(
        n_splits=CALIBRATION_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    return CalibratedClassifierCV(
        estimator=base_model,
        method=CALIBRATION_METHOD,
        cv=calibration_cv,
        n_jobs=1,
    )


# =====================================================
# THRESHOLD SEARCH
# =====================================================

def search_best_threshold(
    y_true,
    probabilities,
) -> tuple[
    float,
    pd.DataFrame,
]:
    """
    Find threshold maximizing balanced accuracy.

    F1 is used as secondary tie-break criterion.
    """

    thresholds = np.arange(
        0.05,
        0.951,
        0.005,
    )

    rows = []

    for threshold in thresholds:

        predictions = (
            probabilities
            >= threshold
        ).astype(int)

        rows.append(
            {
                "threshold":
                    float(threshold),
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
            }
        )

    table = pd.DataFrame(
        rows
    )

    best = (
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
            best["threshold"]
        ),
        table,
    )


# =====================================================
# OUT-OF-FOLD CALIBRATED PREDICTIONS
# =====================================================

def generate_oof_predictions(
    dataset: pd.DataFrame,
    non_methylation_features: list[str],
    methylation_columns: list[str],
) -> pd.DataFrame:
    """
    Generate leakage-aware calibrated OOF probabilities.

    Methylation selection and preprocessing are fitted
    separately within each outer training fold.
    """

    outer_cv = StratifiedKFold(
        n_splits=OUTER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    rows = []

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
            f"FINALIZATION OOF FOLD "
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
            f"Selected methylation features: "
            f"{len(selected_methylation):,}"
        )

        calibrated_model = (
            build_calibrated_model(
                X_train
            )
        )

        print(
            "Training calibrated model..."
        )

        calibrated_model.fit(
            X_train,
            y_train,
        )

        probabilities = (
            calibrated_model
            .predict_proba(
                X_test
            )[:, 1]
        )

        for position in range(
            len(test_df)
        ):

            rows.append(
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
                    "calibrated_probability":
                        float(
                            probabilities[
                                position
                            ]
                        ),
                }
            )

    return pd.DataFrame(
        rows
    )


# =====================================================
# METRICS
# =====================================================

def calculate_final_metrics(
    y_true,
    probabilities,
    threshold: float,
) -> pd.DataFrame:
    """Calculate OOF performance at selected threshold."""

    predictions = (
        probabilities
        >= threshold
    ).astype(int)

    clipped = np.clip(
        probabilities,
        1e-7,
        1 - 1e-7,
    )

    metrics = {
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
    }

    return pd.DataFrame(
        [
            {
                "metric":
                    key,
                "value":
                    value,
            }
            for key, value
            in metrics.items()
        ]
    )


# =====================================================
# FINAL MODEL TRAINING
# =====================================================

def train_final_model(
    dataset: pd.DataFrame,
    non_methylation_features: list[str],
    methylation_columns: list[str],
):
    """
    Train final isotonic-calibrated model
    using all available development data.
    """

    print(
        "\nSelecting final methylation features..."
    )

    selected_methylation = (
        select_methylation_features(
            dataset,
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

    X = dataset[
        feature_columns
    ].copy()

    y = dataset[
        TARGET_COLUMN
    ].copy()

    high_cardinality = (
        identify_high_cardinality_columns(
            X
        )
    )

    if high_cardinality:

        print(
            "\nDropping final high-cardinality "
            "categorical predictors:"
        )

        for column in high_cardinality:
            print(
                f"  - {column}"
            )

        X = X.drop(
            columns=high_cardinality
        )

    final_feature_columns = (
        X.columns.tolist()
    )

    print(
        f"\nFinal raw predictors: "
        f"{len(final_feature_columns):,}"
    )

    print(
        f"Final methylation predictors: "
        f"{len(selected_methylation):,}"
    )

    final_model = (
        build_calibrated_model(
            X
        )
    )

    print(
        "\nTraining final isotonic-calibrated model "
        "on all development patients..."
    )

    final_model.fit(
        X,
        y,
    )

    return (
        final_model,
        final_feature_columns,
        selected_methylation,
        high_cardinality,
    )


# =====================================================
# REPORT
# =====================================================

def create_report(
    dataset: pd.DataFrame,
    metrics: pd.DataFrame,
    threshold: float,
    final_feature_columns: list[str],
    selected_methylation: list[str],
) -> str:
    """Create final model report."""

    metric_map = dict(
        zip(
            metrics["metric"],
            metrics["value"],
        )
    )

    report_lines = [
        "# Final Epigenomic Classification Model",
        "",
        "## Model specification",
        "",
        (
            "- Algorithm: Random Forest classifier"
        ),
        (
            "- Probability calibration: isotonic regression"
        ),
        (
            f"- Random Forest estimators: "
            f"{N_ESTIMATORS}"
        ),
        (
            f"- Calibration folds: "
            f"{CALIBRATION_FOLDS}"
        ),
        "",
        "## Cohort",
        "",
        (
            f"- Development patients: "
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
        "## Predictors",
        "",
        (
            f"- Final raw predictors: "
            f"{len(final_feature_columns):,}"
        ),
        (
            f"- Selected methylation predictors: "
            f"{len(selected_methylation):,}"
        ),
        "",
        "## Development threshold",
        "",
        (
            "- Threshold was selected from leakage-aware "
            "out-of-fold calibrated probabilities."
        ),
        (
            f"- Selected probability threshold: "
            f"{threshold:.3f}"
        ),
        "",
        "## Out-of-fold performance",
        "",
        (
            f"- ROC AUC: "
            f"{metric_map['roc_auc']:.3f}"
        ),
        (
            f"- PR AUC: "
            f"{metric_map['pr_auc']:.3f}"
        ),
        (
            f"- Brier score: "
            f"{metric_map['brier_score']:.4f}"
        ),
        (
            f"- Log loss: "
            f"{metric_map['log_loss']:.4f}"
        ),
        (
            f"- Accuracy: "
            f"{metric_map['accuracy']:.3f}"
        ),
        (
            f"- Balanced accuracy: "
            f"{metric_map['balanced_accuracy']:.3f}"
        ),
        (
            f"- Precision: "
            f"{metric_map['precision']:.3f}"
        ),
        (
            f"- Recall: "
            f"{metric_map['recall']:.3f}"
        ),
        (
            f"- F1-score: "
            f"{metric_map['f1_score']:.3f}"
        ),
        "",
        "## Methodological safeguards",
        "",
        (
            "- Patient/sample identifiers were excluded "
            "from model predictors."
        ),
        (
            "- Survival and direct outcome-proxy variables "
            "were excluded."
        ),
        (
            "- Methylation selection for OOF predictions "
            "was performed within each training fold."
        ),
        (
            "- Isotonic calibration was fitted without "
            "access to outer-test outcomes."
        ),
        "",
        "## Interpretation",
        "",
        (
            "This model is the final development-stage "
            "clinical-genomic-epigenomic classifier."
        ),
        "",
        (
            "The probability threshold is a development-stage "
            "operating threshold derived from internal "
            "cross-validation. It should not be interpreted "
            "as an externally validated clinical decision "
            "threshold."
        ),
        "",
        (
            "Independent external validation is required "
            "before any clinical interpretation or deployment."
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
    """Finalize calibrated epigenomic model."""

    print(
        "Loading Month 1 outcome dataset..."
    )

    baseline = load_dataset(
        BASELINE_FILE
    )

    print(
        "Loading integrated epigenomic dataset..."
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
        f"\nSafe non-methylation predictors: "
        f"{len(non_methylation_features):,}"
    )

    print(
        f"Available methylation predictors: "
        f"{len(methylation_columns):,}"
    )

    if not methylation_columns:
        raise ValueError(
            "No methylation predictors found."
        )

    # =================================================
    # OOF CALIBRATED PROBABILITIES
    # =================================================

    print(
        "\nGenerating leakage-aware "
        "calibrated OOF probabilities..."
    )

    oof_predictions = (
        generate_oof_predictions(
            dataset=dataset,
            non_methylation_features=(
                non_methylation_features
            ),
            methylation_columns=(
                methylation_columns
            ),
        )
    )

    oof_predictions.to_csv(
        OOF_PREDICTIONS_FILE,
        index=False,
    )

    # =================================================
    # THRESHOLD
    # =================================================

    print(
        "\nSearching threshold for "
        "calibrated probabilities..."
    )

    best_threshold, (
        threshold_search
    ) = search_best_threshold(
        y_true=oof_predictions[
            "true_outcome"
        ],
        probabilities=oof_predictions[
            "calibrated_probability"
        ],
    )

    threshold_search.to_csv(
        THRESHOLD_SEARCH_FILE,
        index=False,
    )

    print(
        f"\nFinal calibrated development "
        f"threshold: {best_threshold:.3f}"
    )

    threshold_table = pd.DataFrame(
        {
            "parameter": [
                "calibration_method",
                "threshold_objective",
                "final_threshold",
            ],
            "value": [
                CALIBRATION_METHOD,
                "maximize_balanced_accuracy",
                best_threshold,
            ],
        }
    )

    threshold_table.to_csv(
        FINAL_THRESHOLD_FILE,
        index=False,
    )

    # =================================================
    # METRICS
    # =================================================

    final_metrics = (
        calculate_final_metrics(
            y_true=oof_predictions[
                "true_outcome"
            ],
            probabilities=oof_predictions[
                "calibrated_probability"
            ],
            threshold=best_threshold,
        )
    )

    final_metrics.to_csv(
        FINAL_METRICS_FILE,
        index=False,
    )

    print(
        "\nOOF final metrics:"
    )

    print(
        final_metrics.to_string(
            index=False
        )
    )

    # =================================================
    # TRAIN FINAL MODEL
    # =================================================

    (
        final_model,
        final_feature_columns,
        selected_methylation,
        high_cardinality,
    ) = train_final_model(
        dataset=dataset,
        non_methylation_features=(
            non_methylation_features
        ),
        methylation_columns=(
            methylation_columns
        ),
    )

    # =================================================
    # FEATURE TABLE
    # =================================================

    feature_table = pd.DataFrame(
        {
            "feature":
                final_feature_columns,
        }
    )

    feature_table[
        "feature_type"
    ] = np.where(
        feature_table[
            "feature"
        ].str.startswith(
            "meth_"
        ),
        "methylation",
        "clinical_genomic",
    )

    feature_table.to_csv(
        FINAL_FEATURES_FILE,
        index=False,
    )

    # =================================================
    # SAVE MODEL PACKAGE
    # =================================================

    model_package = {
        "model":
            final_model,
        "patient_column":
            PATIENT_COLUMN,
        "target_column":
            TARGET_COLUMN,
        "calibration_method":
            CALIBRATION_METHOD,
        "threshold":
            best_threshold,
        "feature_columns":
            final_feature_columns,
        "selected_methylation_features":
            selected_methylation,
        "dropped_high_cardinality":
            high_cardinality,
        "random_state":
            RANDOM_STATE,
        "n_estimators":
            N_ESTIMATORS,
        "top_methylation_features":
            TOP_METHYLATION_FEATURES,
    }

    joblib.dump(
        model_package,
        FINAL_MODEL_FILE,
    )

    # =================================================
    # JSON CONFIG
    # =================================================

    metric_dict = {
        row["metric"]:
            float(row["value"])
        for _, row
        in final_metrics.iterrows()
    }

    config = {
        "model_name":
            "Final clinical-genomic-epigenomic Random Forest",
        "calibration":
            CALIBRATION_METHOD,
        "development_threshold":
            float(best_threshold),
        "threshold_objective":
            "maximize_balanced_accuracy",
        "patients":
            int(len(dataset)),
        "positive_outcomes":
            int(
                dataset[
                    TARGET_COLUMN
                ].sum()
            ),
        "negative_outcomes":
            int(
                (
                    dataset[
                        TARGET_COLUMN
                    ] == 0
                ).sum()
            ),
        "raw_predictor_count":
            int(
                len(
                    final_feature_columns
                )
            ),
        "selected_methylation_count":
            int(
                len(
                    selected_methylation
                )
            ),
        "random_forest_estimators":
            N_ESTIMATORS,
        "random_state":
            RANDOM_STATE,
        "internal_oof_metrics":
            metric_dict,
        "external_validation":
            False,
        "clinical_use":
            False,
    }

    FINAL_CONFIG_FILE.write_text(
        json.dumps(
            config,
            indent=2,
        ),
        encoding="utf-8",
    )

    # =================================================
    # REPORT
    # =================================================

    report = create_report(
        dataset=dataset,
        metrics=final_metrics,
        threshold=best_threshold,
        final_feature_columns=(
            final_feature_columns
        ),
        selected_methylation=(
            selected_methylation
        ),
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # =================================================
    # DONE
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "FINAL MODEL COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        f"\nCalibration method: "
        f"{CALIBRATION_METHOD}"
    )

    print(
        f"Final calibrated threshold: "
        f"{best_threshold:.3f}"
    )

    print(
        "\nFinished successfully."
    )

    print(
        "\nFinal model saved to:\n"
        f"{FINAL_MODEL_FILE}"
    )

    print(
        "\nModel configuration saved to:\n"
        f"{FINAL_CONFIG_FILE}"
    )

    print(
        "\nFinal metrics saved to:\n"
        f"{FINAL_METRICS_FILE}"
    )

    print(
        "\nFinal feature list saved to:\n"
        f"{FINAL_FEATURES_FILE}"
    )

    print(
        "\nFinal threshold saved to:\n"
        f"{FINAL_THRESHOLD_FILE}"
    )

    print(
        "\nOOF predictions saved to:\n"
        f"{OOF_PREDICTIONS_FILE}"
    )

    print(
        "\nThreshold search saved to:\n"
        f"{THRESHOLD_SEARCH_FILE}"
    )

    print(
        "\nFinal report saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()