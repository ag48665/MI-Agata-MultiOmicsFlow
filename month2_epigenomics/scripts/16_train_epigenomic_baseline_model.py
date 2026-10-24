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
# FUNCTIONS
# =====================================

def load_dataset(
    file_path: Path,
) -> pd.DataFrame:
    """Load a model feature table."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

    return pd.read_csv(
        file_path,
        low_memory=False,
    )


def validate_dataset(
    df: pd.DataFrame,
    dataset_name: str,
) -> None:
    """Validate required model columns."""

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


def standardize_dataset(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Standardize patient identifiers and target values."""

    df = df.copy()

    df[PATIENT_COLUMN] = (
        df[PATIENT_COLUMN]
        .astype("string")
        .str.strip()
        .str.upper()
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
    """
    Restrict both datasets to the same patients.

    This ensures a fair Month 1 versus Month 2 comparison.
    """

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
            "between the two feature tables."
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

    inconsistent_outcomes = outcome_check[
        outcome_check[
            f"{TARGET_COLUMN}_baseline"
        ]
        != outcome_check[
            f"{TARGET_COLUMN}_epigenomic"
        ]
    ]

    if not inconsistent_outcomes.empty:
        raise ValueError(
            "Outcome values differ between "
            "the Month 1 and Month 2 tables."
        )

    return (
        baseline_aligned,
        epigenomic_aligned,
    )


def create_patient_split(
    epigenomic: pd.DataFrame,
) -> tuple[list[str], list[str]]:
    """Create one shared train/test split."""

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
    """Apply the shared patient split to one dataset."""

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

    numeric_transformer = Pipeline(
        steps=[
            (
                "imputer",
                SimpleImputer(
                    strategy="median"
                ),
            ),
        ]
    )

    categorical_transformer = Pipeline(
        steps=[
            (
                "imputer",
                SimpleImputer(
                    strategy="most_frequent"
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

    preprocessor = ColumnTransformer(
        transformers=[
            (
                "numeric",
                numeric_transformer,
                numeric_features,
            ),
            (
                "categorical",
                categorical_transformer,
                categorical_features,
            ),
        ],
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


def calculate_metrics(
    y_true: pd.Series,
    y_pred: np.ndarray,
    y_probability: np.ndarray,
) -> dict:
    """Calculate classification metrics."""

    metrics = {
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
    }

    if y_true.nunique() == 2:
        metrics["roc_auc"] = (
            roc_auc_score(
                y_true,
                y_probability,
            )
        )
    else:
        metrics["roc_auc"] = np.nan

    return metrics


def save_confusion_matrix(
    y_true: pd.Series,
    y_pred: np.ndarray,
    title: str,
    output_file: Path,
) -> None:
    """Save a confusion matrix figure."""

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
    """Save Random Forest feature importance."""

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
) -> tuple[Pipeline, dict]:
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

    probabilities = (
        pipeline.predict_proba(
            X_test
        )[:, 1]
    )

    metrics = calculate_metrics(
        y_true=y_test,
        y_pred=predictions,
        y_probability=probabilities,
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
        title=f"{model_name}: top feature importances",
    )

    joblib.dump(
        pipeline,
        model_file,
    )

    print(
        f"{model_name} ROC AUC: "
        f"{metrics['roc_auc']:.3f}"
    )

    print(
        f"{model_name} F1-score: "
        f"{metrics['f1_score']:.3f}"
    )

    return pipeline, metrics


def create_comparison_table(
    baseline_metrics: dict,
    epigenomic_metrics: dict,
) -> pd.DataFrame:
    """Compare Month 1 and Month 2 model results."""

    rows = []

    for model_name, metrics in [
        (
            "clinical_genomic",
            baseline_metrics,
        ),
        (
            "clinical_genomic_epigenomic",
            epigenomic_metrics,
        ),
    ]:
        row = {
            "model": model_name,
        }

        row.update(metrics)
        rows.append(row)

    comparison = pd.DataFrame(
        rows
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
    """Create a Markdown model comparison report."""

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
            "- Clinical-genomic-epigenomic input features: "
            f"{epigenomic_feature_count:,}"
        ),
        "",
        "## Model",
        "",
        "- Algorithm: Random Forest Classifier",
        f"- Trees: {N_ESTIMATORS:,}",
        "- Class weighting: balanced",
        (
            "- Both models used the same "
            "patient-level train/test split."
        ),
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
        "## Epigenomic model",
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
            "The model predicts overall survival status. "
            "Because this outcome is imbalanced and does not "
            "include survival time or censoring duration, "
            "the results should be interpreted as a baseline "
            "classification experiment rather than a complete "
            "survival analysis."
        ),
        "",
    ]

    return "\n".join(report_lines)


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Train and compare baseline multi-omics models."""

    print("Loading Month 1 clinical-genomic table...")

    baseline = load_dataset(
        BASELINE_INPUT_FILE
    )

    print("Loading Month 2 epigenomic table...")

    epigenomic = load_dataset(
        EPIGENOMIC_INPUT_FILE
    )

    validate_dataset(
        baseline,
        "Clinical-genomic dataset",
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
        f"Common patients: "
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

    _, baseline_metrics = train_and_evaluate(
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

    _, epigenomic_metrics = train_and_evaluate(
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

    report = create_technical_report(
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
        report,
        encoding="utf-8",
    )

    print("\nModel comparison:")
    print(comparison.to_string(index=False))

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