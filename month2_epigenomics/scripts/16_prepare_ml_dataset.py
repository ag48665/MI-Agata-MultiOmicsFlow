from pathlib import Path

import joblib
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
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
)

MODELS_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "models"
)

FIGURES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "figures"
)

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

for directory in [
    MODELS_DIR,
    FIGURES_DIR,
    TABLES_DIR,
    REPORTS_DIR,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

BASELINE_INPUT_FILE = (
    PROCESSED_DIR
    / "ml_feature_table.csv"
)

EPIGENOMIC_INPUT_FILE = (
    PROCESSED_DIR
    / "ml_epigenomic_feature_table.csv"
)

BASELINE_MODEL_FILE = (
    MODELS_DIR
    / "clinical_genomic_random_forest.pkl"
)

EPIGENOMIC_MODEL_FILE = (
    MODELS_DIR
    / "clinical_genomic_epigenomic_random_forest.pkl"
)

COMPARISON_FILE = (
    TABLES_DIR
    / "epigenomic_model_comparison.csv"
)

BASELINE_REPORT_FILE = (
    TABLES_DIR
    / "clinical_genomic_classification_report.csv"
)

EPIGENOMIC_REPORT_FILE = (
    TABLES_DIR
    / "epigenomic_classification_report.csv"
)

BASELINE_IMPORTANCE_FILE = (
    TABLES_DIR
    / "clinical_genomic_feature_importance.csv"
)

EPIGENOMIC_IMPORTANCE_FILE = (
    TABLES_DIR
    / "epigenomic_feature_importance.csv"
)

TECHNICAL_REPORT_FILE = (
    REPORTS_DIR
    / "epigenomic_baseline_model_report.md"
)


# =====================================
# SETTINGS
# =====================================

PATIENT_COLUMN = "patient_id"
TARGET_COLUMN = "outcome"

TEST_SIZE = 0.20
RANDOM_STATE = 42
N_ESTIMATORS = 500


# =====================================
# DATA FUNCTIONS
# =====================================

def load_dataset(
    file_path: Path,
) -> pd.DataFrame:
    """Load a CSV feature table."""

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
    """Standardize TCGA patient identifiers."""

    if PATIENT_COLUMN not in df.columns:
        raise KeyError(
            f"Column '{PATIENT_COLUMN}' was not found."
        )

    df = df.copy()

    df[PATIENT_COLUMN] = (
        df[PATIENT_COLUMN]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return df


def add_outcome_to_epigenomic_table(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add the Month 1 outcome column to the Month 2 table.

    The merge is performed using standardized patient IDs.
    """

    if TARGET_COLUMN in epigenomic.columns:
        print(
            "Outcome column already exists "
            "in the epigenomic table."
        )
        return epigenomic

    print(
        "Outcome column not found in "
        "the epigenomic table."
    )

    print(
        "Adding outcome from the Month 1 "
        "clinical-genomic table..."
    )

    if TARGET_COLUMN not in baseline.columns:
        raise KeyError(
            f"Month 1 table is missing "
            f"the '{TARGET_COLUMN}' column."
        )

    baseline = standardize_patient_ids(
        baseline
    )

    epigenomic = standardize_patient_ids(
        epigenomic
    )

    outcome_table = baseline[
        [
            PATIENT_COLUMN,
            TARGET_COLUMN,
        ]
    ].copy()

    outcome_table[TARGET_COLUMN] = pd.to_numeric(
        outcome_table[TARGET_COLUMN],
        errors="coerce",
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

    missing_outcomes = int(
        epigenomic[TARGET_COLUMN]
        .isna()
        .sum()
    )

    print(
        "Patients without an outcome "
        f"after matching: {missing_outcomes:,}"
    )

    return epigenomic


def validate_dataset(
    df: pd.DataFrame,
    dataset_name: str,
) -> None:
    """Check required columns and target classes."""

    required_columns = {
        PATIENT_COLUMN,
        TARGET_COLUMN,
    }

    missing_columns = (
        required_columns.difference(
            df.columns
        )
    )

    if missing_columns:
        raise KeyError(
            f"{dataset_name} is missing required "
            f"columns: {sorted(missing_columns)}"
        )

    target = pd.to_numeric(
        df[TARGET_COLUMN],
        errors="coerce",
    ).dropna()

    unique_classes = sorted(
        target.unique().tolist()
    )

    if len(unique_classes) < 2:
        raise ValueError(
            f"{dataset_name} does not contain "
            "at least two target classes."
        )

    unexpected_classes = set(
        unique_classes
    ).difference(
        {0, 1}
    )

    if unexpected_classes:
        raise ValueError(
            f"{dataset_name} contains unexpected "
            f"outcome values: {sorted(unexpected_classes)}"
        )


def standardize_dataset(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Standardize IDs, target values and duplicate records."""

    df = standardize_patient_ids(
        df
    )

    df[TARGET_COLUMN] = pd.to_numeric(
        df[TARGET_COLUMN],
        errors="coerce",
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

    df = (
        df.drop_duplicates(
            subset=[PATIENT_COLUMN],
            keep="first",
        )
        .reset_index(drop=True)
    )

    return df


def align_datasets(
    baseline: pd.DataFrame,
    epigenomic: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict both datasets to the same patients."""

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
            "No common patients were found "
            "between Month 1 and Month 2 tables."
        )

    baseline_aligned = (
        baseline[
            baseline[PATIENT_COLUMN].isin(
                common_patients
            )
        ]
        .sort_values(PATIENT_COLUMN)
        .reset_index(drop=True)
    )

    epigenomic_aligned = (
        epigenomic[
            epigenomic[PATIENT_COLUMN].isin(
                common_patients
            )
        ]
        .sort_values(PATIENT_COLUMN)
        .reset_index(drop=True)
    )

    outcome_check = baseline_aligned[
        [
            PATIENT_COLUMN,
            TARGET_COLUMN,
        ]
    ].merge(
        epigenomic_aligned[
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

    inconsistent = outcome_check[
        outcome_check[
            f"{TARGET_COLUMN}_baseline"
        ]
        != outcome_check[
            f"{TARGET_COLUMN}_epigenomic"
        ]
    ]

    if not inconsistent.empty:
        raise ValueError(
            "Outcome values differ between "
            "the Month 1 and Month 2 tables."
        )

    return (
        baseline_aligned,
        epigenomic_aligned,
    )


# =====================================
# SPLITTING AND PREPROCESSING
# =====================================

def create_patient_split(
    epigenomic: pd.DataFrame,
) -> tuple[list[str], list[str]]:
    """Create one shared stratified train/test split."""

    patient_ids = epigenomic[
        PATIENT_COLUMN
    ]

    target = epigenomic[
        TARGET_COLUMN
    ]

    train_ids, test_ids = train_test_split(
        patient_ids,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=target,
    )

    return (
        train_ids.tolist(),
        test_ids.tolist(),
    )


def split_dataset(
    df: pd.DataFrame,
    train_ids: list[str],
    test_ids: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.Series,
    pd.Series,
]:
    """Apply the same patient split to one dataset."""

    train_df = df[
        df[PATIENT_COLUMN].isin(
            train_ids
        )
    ].copy()

    test_df = df[
        df[PATIENT_COLUMN].isin(
            test_ids
        )
    ].copy()

    feature_columns = [
        column
        for column in df.columns
        if column not in {
            PATIENT_COLUMN,
            TARGET_COLUMN,
        }
    ]

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

    return (
        X_train,
        X_test,
        y_train,
        y_test,
    )


def build_pipeline(
    X_train: pd.DataFrame,
) -> Pipeline:
    """Create preprocessing and Random Forest pipeline."""

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

    print(
        f"Numeric features: "
        f"{len(numeric_features):,}"
    )

    print(
        f"Categorical features: "
        f"{len(categorical_features):,}"
    )

    transformers = []

    if numeric_features:
        numeric_transformer = Pipeline(
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
                numeric_transformer,
                numeric_features,
            )
        )

    if categorical_features:
        categorical_transformer = Pipeline(
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
                categorical_transformer,
                categorical_features,
            )
        )

    if not transformers:
        raise ValueError(
            "No model features were detected."
        )

    preprocessor = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
    )

    classifier = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        class_weight="balanced",
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


# =====================================
# EVALUATION
# =====================================

def calculate_metrics(
    y_true: pd.Series,
    y_pred: np.ndarray,
    y_probability: np.ndarray,
) -> dict[str, float]:
    """Calculate classification metrics."""

    return {
        "accuracy": accuracy_score(
            y_true,
            y_pred,
        ),
        "balanced_accuracy": (
            balanced_accuracy_score(
                y_true,
                y_pred,
            )
        ),
        "precision": precision_score(
            y_true,
            y_pred,
            zero_division=0,
        ),
        "recall": recall_score(
            y_true,
            y_pred,
            zero_division=0,
        ),
        "f1_score": f1_score(
            y_true,
            y_pred,
            zero_division=0,
        ),
        "roc_auc": roc_auc_score(
            y_true,
            y_probability,
        ),
    }


def save_confusion_matrix(
    y_true: pd.Series,
    y_pred: np.ndarray,
    title: str,
    output_file: Path,
) -> None:
    """Save confusion matrix figure."""

    matrix = confusion_matrix(
        y_true,
        y_pred,
    )

    display = ConfusionMatrixDisplay(
        confusion_matrix=matrix,
    )

    display.plot()

    plt.title(title)
    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def save_feature_importance(
    pipeline: Pipeline,
    output_table: Path,
    output_figure: Path,
    title: str,
) -> pd.DataFrame:
    """Save model feature-importance results."""

    feature_names = (
        pipeline
        .named_steps["preprocessor"]
        .get_feature_names_out()
    )

    importances = (
        pipeline
        .named_steps["classifier"]
        .feature_importances_
    )

    if len(feature_names) != len(importances):
        raise ValueError(
            "Feature-name count does not match "
            "the feature-importance count."
        )

    importance_table = (
        pd.DataFrame(
            {
                "feature": feature_names,
                "importance": importances,
            }
        )
        .sort_values(
            "importance",
            ascending=False,
        )
        .reset_index(drop=True)
    )

    importance_table[
        "rank"
    ] = importance_table.index + 1

    importance_table.to_csv(
        output_table,
        index=False,
    )

    top_features = (
        importance_table
        .head(20)
        .sort_values(
            "importance",
            ascending=True,
        )
    )

    plt.figure(
        figsize=(10, 8)
    )

    plt.barh(
        top_features["feature"],
        top_features["importance"],
    )

    plt.xlabel("Feature importance")
    plt.title(title)
    plt.tight_layout()

    plt.savefig(
        output_figure,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    return importance_table


def train_and_evaluate(
    model_name: str,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: pd.Series,
    y_test: pd.Series,
    model_file: Path,
    classification_report_file: Path,
    importance_file: Path,
    confusion_figure: Path,
    importance_figure: Path,
) -> dict[str, float]:
    """Train and evaluate one Random Forest model."""

    print(
        f"\nTraining {model_name}..."
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

    probabilities = pipeline.predict_proba(
        X_test
    )

    classes = (
        pipeline
        .named_steps["classifier"]
        .classes_
    )

    positive_class_index = int(
        np.where(
            classes == 1
        )[0][0]
    )

    positive_probabilities = probabilities[
        :,
        positive_class_index,
    ]

    metrics = calculate_metrics(
        y_true=y_test,
        y_pred=predictions,
        y_probability=positive_probabilities,
    )

    report = classification_report(
        y_test,
        predictions,
        output_dict=True,
        zero_division=0,
    )

    pd.DataFrame(
        report
    ).transpose().to_csv(
        classification_report_file
    )

    save_confusion_matrix(
        y_true=y_test,
        y_pred=predictions,
        title=f"{model_name} confusion matrix",
        output_file=confusion_figure,
    )

    save_feature_importance(
        pipeline=pipeline,
        output_table=importance_file,
        output_figure=importance_figure,
        title=f"{model_name}: feature importance",
    )

    joblib.dump(
        pipeline,
        model_file,
    )

    print(
        f"{model_name} accuracy: "
        f"{metrics['accuracy']:.3f}"
    )

    print(
        f"{model_name} balanced accuracy: "
        f"{metrics['balanced_accuracy']:.3f}"
    )

    print(
        f"{model_name} F1-score: "
        f"{metrics['f1_score']:.3f}"
    )

    print(
        f"{model_name} ROC AUC: "
        f"{metrics['roc_auc']:.3f}"
    )

    return metrics


# =====================================
# COMPARISON AND REPORT
# =====================================

def create_comparison_table(
    baseline_metrics: dict[str, float],
    epigenomic_metrics: dict[str, float],
) -> pd.DataFrame:
    """Create model comparison table."""

    comparison = pd.DataFrame(
        [
            {
                "model": "clinical_genomic",
                **baseline_metrics,
            },
            {
                "model": (
                    "clinical_genomic_epigenomic"
                ),
                **epigenomic_metrics,
            },
        ]
    )

    baseline_roc = comparison.loc[
        comparison["model"]
        == "clinical_genomic",
        "roc_auc",
    ].iloc[0]

    comparison[
        "roc_auc_change_vs_baseline"
    ] = (
        comparison["roc_auc"]
        - baseline_roc
    )

    return comparison


def create_technical_report(
    comparison: pd.DataFrame,
    patient_count: int,
    train_count: int,
    test_count: int,
    baseline_feature_count: int,
    epigenomic_feature_count: int,
) -> str:
    """Create Markdown report for both models."""

    baseline_row = comparison[
        comparison["model"]
        == "clinical_genomic"
    ].iloc[0]

    epigenomic_row = comparison[
        comparison["model"]
        == "clinical_genomic_epigenomic"
    ].iloc[0]

    roc_change = (
        epigenomic_row["roc_auc"]
        - baseline_row["roc_auc"]
    )

    report_lines = [
        "# Epigenomic Baseline Model Report",
        "",
        "## Cohort",
        "",
        f"- Matched patients: {patient_count:,}",
        f"- Training patients: {train_count:,}",
        f"- Testing patients: {test_count:,}",
        (
            "- Clinical-genomic input features: "
            f"{baseline_feature_count:,}"
        ),
        (
            "- Epigenomic input features: "
            f"{epigenomic_feature_count:,}"
        ),
        "",
        "## Model settings",
        "",
        "- Algorithm: Random Forest Classifier",
        f"- Trees: {N_ESTIMATORS:,}",
        "- Class weight: balanced",
        "- Train/test split: 80/20",
        "- Shared patient split for both models",
        "",
        "## Clinical-genomic model",
        "",
        f"- Accuracy: {baseline_row['accuracy']:.3f}",
        (
            "- Balanced accuracy: "
            f"{baseline_row['balanced_accuracy']:.3f}"
        ),
        f"- Precision: {baseline_row['precision']:.3f}",
        f"- Recall: {baseline_row['recall']:.3f}",
        f"- F1-score: {baseline_row['f1_score']:.3f}",
        f"- ROC AUC: {baseline_row['roc_auc']:.3f}",
        "",
        "## Clinical-genomic-epigenomic model",
        "",
        f"- Accuracy: {epigenomic_row['accuracy']:.3f}",
        (
            "- Balanced accuracy: "
            f"{epigenomic_row['balanced_accuracy']:.3f}"
        ),
        f"- Precision: {epigenomic_row['precision']:.3f}",
        f"- Recall: {epigenomic_row['recall']:.3f}",
        f"- F1-score: {epigenomic_row['f1_score']:.3f}",
        f"- ROC AUC: {epigenomic_row['roc_auc']:.3f}",
        "",
        "## Epigenomic contribution",
        "",
        (
            "- ROC AUC change after adding methylation: "
            f"{roc_change:+.3f}"
        ),
        "",
        "## Interpretation note",
        "",
        (
            "The target is overall survival status. "
            "This is a baseline classification experiment "
            "and not a full time-to-event survival analysis."
        ),
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Train and compare Month 1 and Month 2 models."""

    print(
        "Loading Month 1 clinical-genomic table..."
    )

    baseline = load_dataset(
        BASELINE_INPUT_FILE
    )

    print(
        "Loading Month 2 epigenomic table..."
    )

    epigenomic = load_dataset(
        EPIGENOMIC_INPUT_FILE
    )

    baseline = standardize_patient_ids(
        baseline
    )

    epigenomic = standardize_patient_ids(
        epigenomic
    )

    validate_dataset(
        baseline,
        "Clinical-genomic dataset",
    )

    epigenomic = add_outcome_to_epigenomic_table(
        baseline=baseline,
        epigenomic=epigenomic,
    )

    validate_dataset(
        epigenomic,
        "Epigenomic dataset",
    )

    baseline = standardize_dataset(
        baseline
    )

    epigenomic = standardize_dataset(
        epigenomic
    )

    baseline, epigenomic = align_datasets(
        baseline,
        epigenomic,
    )

    print(
        f"\nCommon patients: "
        f"{len(epigenomic):,}"
    )

    print("\nOutcome distribution:")

    print(
        epigenomic[
            TARGET_COLUMN
        ].value_counts(
            dropna=False
        )
    )

    train_ids, test_ids = (
        create_patient_split(
            epigenomic
        )
    )

    (
        X_train_baseline,
        X_test_baseline,
        y_train_baseline,
        y_test_baseline,
    ) = split_dataset(
        baseline,
        train_ids,
        test_ids,
    )

    (
        X_train_epigenomic,
        X_test_epigenomic,
        y_train_epigenomic,
        y_test_epigenomic,
    ) = split_dataset(
        epigenomic,
        train_ids,
        test_ids,
    )

    print(
        f"\nTraining patients: "
        f"{len(y_train_epigenomic):,}"
    )

    print(
        f"Testing patients: "
        f"{len(y_test_epigenomic):,}"
    )

    baseline_metrics = train_and_evaluate(
        model_name="Clinical-genomic model",
        X_train=X_train_baseline,
        X_test=X_test_baseline,
        y_train=y_train_baseline,
        y_test=y_test_baseline,
        model_file=BASELINE_MODEL_FILE,
        classification_report_file=(
            BASELINE_REPORT_FILE
        ),
        importance_file=(
            BASELINE_IMPORTANCE_FILE
        ),
        confusion_figure=(
            FIGURES_DIR
            / "clinical_genomic_confusion_matrix.png"
        ),
        importance_figure=(
            FIGURES_DIR
            / "clinical_genomic_feature_importance.png"
        ),
    )

    epigenomic_metrics = train_and_evaluate(
        model_name=(
            "Clinical-genomic-epigenomic model"
        ),
        X_train=X_train_epigenomic,
        X_test=X_test_epigenomic,
        y_train=y_train_epigenomic,
        y_test=y_test_epigenomic,
        model_file=EPIGENOMIC_MODEL_FILE,
        classification_report_file=(
            EPIGENOMIC_REPORT_FILE
        ),
        importance_file=(
            EPIGENOMIC_IMPORTANCE_FILE
        ),
        confusion_figure=(
            FIGURES_DIR
            / "epigenomic_confusion_matrix.png"
        ),
        importance_figure=(
            FIGURES_DIR
            / "epigenomic_feature_importance.png"
        ),
    )

    comparison = create_comparison_table(
        baseline_metrics,
        epigenomic_metrics,
    )

    comparison.to_csv(
        COMPARISON_FILE,
        index=False,
    )

    technical_report = create_technical_report(
        comparison=comparison,
        patient_count=len(epigenomic),
        train_count=len(y_train_epigenomic),
        test_count=len(y_test_epigenomic),
        baseline_feature_count=(
            X_train_baseline.shape[1]
        ),
        epigenomic_feature_count=(
            X_train_epigenomic.shape[1]
        ),
    )

    TECHNICAL_REPORT_FILE.write_text(
        technical_report,
        encoding="utf-8",
    )

    print("\nModel comparison:")
    print(
        comparison.to_string(
            index=False
        )
    )

    print("\nFinished successfully.")

    print(
        "Model comparison saved to:\n"
        f"{COMPARISON_FILE}"
    )

    print(
        "\nTechnical report saved to:\n"
        f"{TECHNICAL_REPORT_FILE}"
    )


if __name__ == "__main__":
    main()