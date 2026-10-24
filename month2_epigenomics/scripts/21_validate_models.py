from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pandas.api.types import is_numeric_dtype

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
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
    / "model_validation_folds.csv"
)

SUMMARY_FILE = (
    TABLES_DIR
    / "model_validation_summary.csv"
)

LEAKAGE_AUDIT_FILE = (
    TABLES_DIR
    / "model_leakage_audit.csv"
)

SELECTED_FEATURES_FILE = (
    TABLES_DIR
    / "cv_selected_methylation_features.csv"
)

FIGURE_FILE = (
    FIGURES_DIR
    / "model_validation_roc_auc.png"
)

REPORT_FILE = (
    REPORTS_DIR
    / "model_validation_report.md"
)


# =====================================================
# SETTINGS
# =====================================================

PATIENT_COLUMN = "patient_id"
TARGET_COLUMN = "outcome"

N_SPLITS = 5
RANDOM_STATE = 42

TOP_METHYLATION_FEATURES = 500

N_ESTIMATORS = 300

MAX_CATEGORICAL_LEVELS = 30

MIN_VARIANCE = 1e-8


# =====================================================
# LOADING
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
# TARGET
# =====================================================

def add_outcome_to_epigenomic(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
) -> pd.DataFrame:
    """Add outcome to epigenomic data if missing."""

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
            f"Month 1 table is missing "
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

    missing = int(
        epigenomic[
            TARGET_COLUMN
        ]
        .isna()
        .sum()
    )

    print(
        f"Missing outcomes after matching: "
        f"{missing:,}"
    )

    return epigenomic


def clean_target(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Clean binary classification target."""

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
# ALIGN PATIENTS
# =====================================================

def align_datasets(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Restrict both datasets to identical patients."""

    common_patients = sorted(
        set(
            baseline[PATIENT_COLUMN]
        ).intersection(
            set(
                epigenomic[PATIENT_COLUMN]
            )
        )
    )

    if not common_patients:
        raise ValueError(
            "No common patients found."
        )

    baseline = (
        baseline[
            baseline[PATIENT_COLUMN].isin(
                common_patients
            )
        ]
        .sort_values(PATIENT_COLUMN)
        .reset_index(drop=True)
    )

    epigenomic = (
        epigenomic[
            epigenomic[PATIENT_COLUMN].isin(
                common_patients
            )
        ]
        .sort_values(PATIENT_COLUMN)
        .reset_index(drop=True)
    )

    comparison = baseline[
        [
            PATIENT_COLUMN,
            TARGET_COLUMN,
        ]
    ].merge(
        epigenomic[
            [
                PATIENT_COLUMN,
                TARGET_COLUMN,
            ]
        ],
        on=PATIENT_COLUMN,
        suffixes=(
            "_baseline",
            "_epigenomic",
        ),
        validate="one_to_one",
    )

    mismatches = comparison[
        comparison[
            f"{TARGET_COLUMN}_baseline"
        ]
        != comparison[
            f"{TARGET_COLUMN}_epigenomic"
        ]
    ]

    if not mismatches.empty:
        raise ValueError(
            "Outcome mismatch between datasets."
        )

    return (
        baseline,
        epigenomic,
    )


# =====================================================
# LEAKAGE FILTER
# =====================================================

def is_forbidden_feature(
    column: str,
) -> tuple[bool, str]:
    """Identify potential leakage or identifier features."""

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
        return True, "exact_exclusion"

    identifier_fragments = [
        "sample_id",
        "patient_id",
        "case_id",
        "barcode",
        "submitter_id",
        "project_id",
        "study_id",
        "tissue_source_site",
    ]

    for fragment in identifier_fragments:
        if fragment in lower:
            return True, (
                f"identifier:{fragment}"
            )

    # Survival leakage
    if lower.endswith("_months"):
        return True, "survival_time"

    if "survival" in lower:
        return True, "survival_variable"

    if lower.startswith("os_"):
        return True, "os_variable"

    if lower.startswith("dfs_"):
        return True, "dfs_variable"

    if lower.startswith("pfs_"):
        return True, "pfs_variable"

    if lower.startswith("dss_"):
        return True, "dss_variable"

    # Direct target-like names
    target_fragments = [
        "death",
        "deceased",
        "alive",
        "living",
        "mortality",
    ]

    for fragment in target_fragments:
        if fragment in lower:
            return True, (
                f"target_proxy:{fragment}"
            )

    return False, ""


def audit_features(
    df: pd.DataFrame,
    dataset_name: str,
) -> tuple[
    list[str],
    pd.DataFrame,
]:
    """Create feature leakage audit."""

    allowed = []

    rows = []

    for column in df.columns:

        forbidden, reason = (
            is_forbidden_feature(
                column
            )
        )

        is_methylation = (
            str(column).startswith(
                "meth_"
            )
        )

        rows.append(
            {
                "dataset":
                    dataset_name,
                "feature":
                    column,
                "is_methylation":
                    is_methylation,
                "excluded":
                    forbidden,
                "reason":
                    reason,
            }
        )

        if not forbidden:
            allowed.append(
                column
            )

    audit = pd.DataFrame(
        rows
    )

    return (
        allowed,
        audit,
    )


# =====================================================
# METHYLATION SELECTION
# =====================================================

def identify_methylation_columns(
    df: pd.DataFrame,
) -> list[str]:
    """Identify methylation feature columns."""

    return [
        column
        for column in df.columns
        if str(column).startswith(
            "meth_"
        )
    ]


def select_fold_methylation_features(
    train_df: pd.DataFrame,
    methylation_columns: list[str],
) -> list[str]:
    """
    Select top methylation features using
    training fold only.
    """

    methylation = train_df[
        methylation_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    variances = methylation.var(
        axis=0,
        skipna=True,
    )

    variances = variances[
        variances.notna()
        & (variances > MIN_VARIANCE)
    ]

    selected = (
        variances
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
# PREPROCESSING
# =====================================================

def remove_high_cardinality_columns(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Remove categorical variables with many levels."""

    train_df = train_df.copy()
    test_df = test_df.copy()

    categorical = [
        column
        for column in train_df.columns
        if not is_numeric_dtype(
            train_df[column]
        )
    ]

    high_cardinality = [
        column
        for column in categorical
        if train_df[column].nunique(
            dropna=True
        ) > MAX_CATEGORICAL_LEVELS
    ]

    if high_cardinality:

        train_df = train_df.drop(
            columns=high_cardinality
        )

        test_df = test_df.drop(
            columns=high_cardinality
        )

    return (
        train_df,
        test_df,
    )


def build_pipeline(
    X_train: pd.DataFrame,
) -> Pipeline:
    """Create leakage-safe preprocessing and RF pipeline."""

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


# =====================================================
# METRICS
# =====================================================

def calculate_metrics(
    y_true: pd.Series,
    predictions: np.ndarray,
    probabilities: np.ndarray,
) -> dict[str, float]:
    """Calculate validation metrics."""

    return {
        "accuracy": accuracy_score(
            y_true,
            predictions,
        ),
        "balanced_accuracy": (
            balanced_accuracy_score(
                y_true,
                predictions,
            )
        ),
        "precision": precision_score(
            y_true,
            predictions,
            zero_division=0,
        ),
        "recall": recall_score(
            y_true,
            predictions,
            zero_division=0,
        ),
        "f1_score": f1_score(
            y_true,
            predictions,
            zero_division=0,
        ),
        "roc_auc": roc_auc_score(
            y_true,
            probabilities,
        ),
    }


# =====================================================
# TRAIN ONE FOLD
# =====================================================

def evaluate_fold(
    model_name: str,
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_columns: list[str],
    fold_number: int,
) -> dict:
    """Train and evaluate one cross-validation fold."""

    X_train = train_df[
        feature_columns
    ].copy()

    X_test = test_df[
        feature_columns
    ].copy()

    y_train = train_df[
        TARGET_COLUMN
    ]

    y_test = test_df[
        TARGET_COLUMN
    ]

    (
        X_train,
        X_test,
    ) = remove_high_cardinality_columns(
        X_train,
        X_test,
    )

    pipeline = build_pipeline(
        X_train
    )

    pipeline.fit(
        X_train,
        y_train,
    )

    predictions = pipeline.predict(
        X_test
    )

    probability_matrix = (
        pipeline.predict_proba(
            X_test
        )
    )

    classes = (
        pipeline
        .named_steps["classifier"]
        .classes_
    )

    positive_index = int(
        np.where(
            classes == 1
        )[0][0]
    )

    probabilities = (
        probability_matrix[
            :,
            positive_index,
        ]
    )

    metrics = calculate_metrics(
        y_true=y_test,
        predictions=predictions,
        probabilities=probabilities,
    )

    metrics.update(
        {
            "model":
                model_name,
            "fold":
                fold_number,
            "training_patients":
                len(train_df),
            "testing_patients":
                len(test_df),
            "raw_features":
                len(feature_columns),
            "features_after_cardinality_qc":
                X_train.shape[1],
            "training_events":
                int(y_train.sum()),
            "testing_events":
                int(y_test.sum()),
        }
    )

    return metrics


# =====================================================
# VALIDATION
# =====================================================

def run_cross_validation(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
    baseline_features: list[str],
    epigenomic_non_methylation: list[str],
    methylation_columns: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Run paired stratified cross-validation."""

    splitter = StratifiedKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    target = baseline[
        TARGET_COLUMN
    ].to_numpy()

    fold_results = []

    selected_feature_rows = []

    for fold_number, (
        train_indices,
        test_indices,
    ) in enumerate(
        splitter.split(
            baseline,
            target,
        ),
        start=1,
    ):

        print(
            f"\n========== Fold "
            f"{fold_number}/{N_SPLITS} =========="
        )

        baseline_train = (
            baseline
            .iloc[
                train_indices
            ]
            .reset_index(drop=True)
        )

        baseline_test = (
            baseline
            .iloc[
                test_indices
            ]
            .reset_index(drop=True)
        )

        epigenomic_train = (
            epigenomic
            .iloc[
                train_indices
            ]
            .reset_index(drop=True)
        )

        epigenomic_test = (
            epigenomic
            .iloc[
                test_indices
            ]
            .reset_index(drop=True)
        )

        # ---------------------------------
        # BASELINE
        # ---------------------------------

        baseline_result = evaluate_fold(
            model_name="clinical_genomic",
            train_df=baseline_train,
            test_df=baseline_test,
            feature_columns=baseline_features,
            fold_number=fold_number,
        )

        fold_results.append(
            baseline_result
        )

        print(
            "Clinical-genomic ROC AUC: "
            f"{baseline_result['roc_auc']:.3f}"
        )

        # ---------------------------------
        # METHYLATION SELECTION
        # TRAINING FOLD ONLY
        # ---------------------------------

        selected_methylation = (
            select_fold_methylation_features(
                train_df=epigenomic_train,
                methylation_columns=(
                    methylation_columns
                ),
            )
        )

        for rank, feature in enumerate(
            selected_methylation,
            start=1,
        ):

            selected_feature_rows.append(
                {
                    "fold":
                        fold_number,
                    "rank":
                        rank,
                    "feature":
                        feature,
                }
            )

        epigenomic_features = (
            epigenomic_non_methylation
            + selected_methylation
        )

        epigenomic_features = list(
            dict.fromkeys(
                epigenomic_features
            )
        )

        epigenomic_result = evaluate_fold(
            model_name=(
                "clinical_genomic_epigenomic"
            ),
            train_df=epigenomic_train,
            test_df=epigenomic_test,
            feature_columns=(
                epigenomic_features
            ),
            fold_number=fold_number,
        )

        fold_results.append(
            epigenomic_result
        )

        print(
            "Epigenomic ROC AUC: "
            f"{epigenomic_result['roc_auc']:.3f}"
        )

    return (
        pd.DataFrame(
            fold_results
        ),
        pd.DataFrame(
            selected_feature_rows
        ),
    )


# =====================================================
# SUMMARY
# =====================================================

def create_summary(
    fold_results: pd.DataFrame,
) -> pd.DataFrame:
    """Create mean/std CV performance summary."""

    metrics = [
        "accuracy",
        "balanced_accuracy",
        "precision",
        "recall",
        "f1_score",
        "roc_auc",
    ]

    rows = []

    for model_name in (
        fold_results[
            "model"
        ].unique()
    ):

        model_data = fold_results[
            fold_results[
                "model"
            ] == model_name
        ]

        row = {
            "model":
                model_name,
            "folds":
                len(model_data),
        }

        for metric in metrics:

            row[
                f"{metric}_mean"
            ] = (
                model_data[
                    metric
                ].mean()
            )

            row[
                f"{metric}_std"
            ] = (
                model_data[
                    metric
                ].std()
            )

        rows.append(
            row
        )

    summary = pd.DataFrame(
        rows
    )

    baseline_auc = summary.loc[
        summary["model"]
        == "clinical_genomic",
        "roc_auc_mean",
    ]

    if not baseline_auc.empty:

        baseline_value = (
            baseline_auc.iloc[0]
        )

        summary[
            "roc_auc_change_vs_baseline"
        ] = (
            summary[
                "roc_auc_mean"
            ]
            - baseline_value
        )

    return summary


# =====================================================
# FIGURE
# =====================================================

def plot_validation_results(
    fold_results: pd.DataFrame,
) -> None:
    """Create fold-level ROC AUC comparison plot."""

    baseline = fold_results[
        fold_results["model"]
        == "clinical_genomic"
    ]

    epigenomic = fold_results[
        fold_results["model"]
        == "clinical_genomic_epigenomic"
    ]

    folds = baseline[
        "fold"
    ].to_numpy()

    plt.figure(
        figsize=(
            9,
            6,
        )
    )

    plt.plot(
        folds,
        baseline[
            "roc_auc"
        ],
        marker="o",
        label=(
            "Clinical + genomic"
        ),
    )

    plt.plot(
        folds,
        epigenomic[
            "roc_auc"
        ],
        marker="o",
        label=(
            "Clinical + genomic + methylation"
        ),
    )

    plt.axhline(
        0.5,
        linestyle="--",
        label="Random discrimination",
    )

    plt.xlabel(
        "Cross-validation fold"
    )

    plt.ylabel(
        "ROC AUC"
    )

    plt.title(
        "Leakage-aware 5-fold model validation"
    )

    plt.ylim(
        0.0,
        1.05,
    )

    plt.xticks(
        folds
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        FIGURE_FILE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# =====================================================
# REPORT
# =====================================================

def dataframe_markdown(
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
    baseline: pd.DataFrame,
    methylation_columns: list[str],
    baseline_features: list[str],
    epigenomic_non_methylation: list[str],
    audit: pd.DataFrame,
    fold_results: pd.DataFrame,
    summary: pd.DataFrame,
    selected_features: pd.DataFrame,
) -> str:
    """Create final validation report."""

    excluded = audit[
        audit[
            "excluded"
        ]
    ]

    methylation_stability = (
        selected_features[
            "feature"
        ]
        .value_counts()
        .rename_axis(
            "feature"
        )
        .reset_index(
            name="folds_selected"
        )
    )

    methylation_stability[
        "selection_frequency"
    ] = (
        methylation_stability[
            "folds_selected"
        ]
        / N_SPLITS
    )

    methylation_stability = (
        methylation_stability
        .head(20)
    )

    report_lines = [
        "# Model Validation and Data Leakage Audit",
        "",
        "## Objective",
        "",
        (
            "This validation stage was created to investigate "
            "the unusually high performance observed in the "
            "initial Month 2 epigenomic classifier."
        ),
        "",
        (
            "The validation uses paired stratified cross-validation "
            "and performs methylation feature selection independently "
            "inside each training fold."
        ),
        "",
        "## Cohort",
        "",
        (
            f"- Matched patients: "
            f"{len(baseline):,}"
        ),
        (
            f"- Outcome-positive patients: "
            f"{int(baseline[TARGET_COLUMN].sum()):,}"
        ),
        (
            f"- Outcome-negative patients: "
            f"{int((baseline[TARGET_COLUMN] == 0).sum()):,}"
        ),
        "",
        "## Cross-validation design",
        "",
        (
            f"- Stratified folds: "
            f"{N_SPLITS}"
        ),
        (
            f"- Random Forest trees per fold: "
            f"{N_ESTIMATORS}"
        ),
        (
            "- Identical patient folds were used for "
            "clinical-genomic and epigenomic models."
        ),
        (
            "- Missing-value imputation was fitted using "
            "training data only."
        ),
        (
            "- Categorical encoding was fitted using "
            "training data only."
        ),
        (
            "- High-cardinality categorical variables "
            "were excluded."
        ),
        (
            "- Methylation variance selection was performed "
            "using training data only."
        ),
        "",
        "## Feature sets",
        "",
        (
            "- Clinical-genomic candidate features: "
            f"{len(baseline_features):,}"
        ),
        (
            "- Non-methylation Month 2 candidate features: "
            f"{len(epigenomic_non_methylation):,}"
        ),
        (
            "- Available methylation features: "
            f"{len(methylation_columns):,}"
        ),
        (
            "- Methylation features selected per fold: "
            f"{TOP_METHYLATION_FEATURES:,}"
        ),
        "",
        "## Leakage audit",
        "",
        (
            f"- Features excluded as identifiers, survival "
            f"variables or target proxies: {len(excluded):,}"
        ),
        "",
        dataframe_markdown(
            excluded[
                [
                    "dataset",
                    "feature",
                    "reason",
                ]
            ].head(30)
        ),
        "",
        "## Cross-validation summary",
        "",
        dataframe_markdown(
            summary.round(4)
        ),
        "",
        "## Fold-level results",
        "",
        dataframe_markdown(
            fold_results[
                [
                    "model",
                    "fold",
                    "accuracy",
                    "balanced_accuracy",
                    "f1_score",
                    "roc_auc",
                ]
            ].round(4)
        ),
        "",
        "## Methylation feature-selection stability",
        "",
        (
            "The following table shows features most frequently "
            "selected among the five independent training folds."
        ),
        "",
        dataframe_markdown(
            methylation_stability
        ),
        "",
        "## Interpretation",
        "",
        (
            "If the cross-validated epigenomic ROC AUC remains "
            "close to 1.0 across independent folds after the "
            "leakage controls above, the strong predictive signal "
            "is more credible."
        ),
        "",
        (
            "If performance decreases substantially relative to "
            "the original hold-out result, the original ROC AUC "
            "of 1.0 should not be reported as a stable estimate "
            "of generalisation performance."
        ),
        "",
        (
            "Balanced accuracy, recall and F1-score should also "
            "be considered because the outcome distribution is "
            "imbalanced."
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
    """Run leakage-aware model validation."""

    print(
        "Loading Month 1 clinical-genomic dataset..."
    )

    baseline = load_dataset(
        BASELINE_FILE
    )

    print(
        "Loading integrated Month 2 epigenomic dataset..."
    )

    epigenomic = load_dataset(
        EPIGENOMIC_FILE
    )

    baseline = clean_target(
        baseline
    )

    epigenomic = add_outcome_to_epigenomic(
        baseline=baseline,
        epigenomic=epigenomic,
    )

    epigenomic = clean_target(
        epigenomic
    )

    (
        baseline,
        epigenomic,
    ) = align_datasets(
        baseline,
        epigenomic,
    )

    print(
        f"\nMatched patients: "
        f"{len(baseline):,}"
    )

    print(
        "\nOutcome distribution:"
    )

    print(
        baseline[
            TARGET_COLUMN
        ]
        .value_counts()
        .sort_index()
    )

    # ---------------------------------
    # LEAKAGE AUDIT
    # ---------------------------------

    (
        baseline_allowed,
        baseline_audit,
    ) = audit_features(
        baseline,
        "clinical_genomic",
    )

    (
        epigenomic_allowed,
        epigenomic_audit,
    ) = audit_features(
        epigenomic,
        "clinical_genomic_epigenomic",
    )

    audit = pd.concat(
        [
            baseline_audit,
            epigenomic_audit,
        ],
        ignore_index=True,
    )

    audit.to_csv(
        LEAKAGE_AUDIT_FILE,
        index=False,
    )

    baseline_features = [
        column
        for column in baseline_allowed
        if column not in {
            PATIENT_COLUMN,
            TARGET_COLUMN,
        }
    ]

    methylation_columns = (
        identify_methylation_columns(
            epigenomic
        )
    )

    epigenomic_non_methylation = [
        column
        for column in epigenomic_allowed
        if (
            column
            not in {
                PATIENT_COLUMN,
                TARGET_COLUMN,
            }
            and not str(
                column
            ).startswith(
                "meth_"
            )
        )
    ]

    print(
        f"\nClinical-genomic predictors "
        f"after leakage audit: "
        f"{len(baseline_features):,}"
    )

    print(
        f"Available methylation features: "
        f"{len(methylation_columns):,}"
    )

    print(
        f"Month 2 non-methylation predictors: "
        f"{len(epigenomic_non_methylation):,}"
    )

    # ---------------------------------
    # CV
    # ---------------------------------

    print(
        "\nStarting leakage-aware "
        "cross-validation..."
    )

    (
        fold_results,
        selected_features,
    ) = run_cross_validation(
        baseline=baseline,
        epigenomic=epigenomic,
        baseline_features=baseline_features,
        epigenomic_non_methylation=(
            epigenomic_non_methylation
        ),
        methylation_columns=(
            methylation_columns
        ),
    )

    # ---------------------------------
    # SAVE
    # ---------------------------------

    fold_results.to_csv(
        FOLD_RESULTS_FILE,
        index=False,
    )

    selected_features.to_csv(
        SELECTED_FEATURES_FILE,
        index=False,
    )

    summary = create_summary(
        fold_results
    )

    summary.to_csv(
        SUMMARY_FILE,
        index=False,
    )

    plot_validation_results(
        fold_results
    )

    report = create_report(
        baseline=baseline,
        methylation_columns=(
            methylation_columns
        ),
        baseline_features=(
            baseline_features
        ),
        epigenomic_non_methylation=(
            epigenomic_non_methylation
        ),
        audit=audit,
        fold_results=fold_results,
        summary=summary,
        selected_features=(
            selected_features
        ),
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # ---------------------------------
    # RESULTS
    # ---------------------------------

    print(
        "\n========================================"
    )

    print(
        "VALIDATION SUMMARY"
    )

    print(
        "========================================"
    )

    print(
        summary.to_string(
            index=False
        )
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
        "\nLeakage audit saved to:\n"
        f"{LEAKAGE_AUDIT_FILE}"
    )

    print(
        "\nSelected methylation features saved to:\n"
        f"{SELECTED_FEATURES_FILE}"
    )

    print(
        "\nValidation figure saved to:\n"
        f"{FIGURE_FILE}"
    )

    print(
        "\nValidation report saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()