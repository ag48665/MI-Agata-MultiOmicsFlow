from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split


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


INPUT_FILE = (
    PROCESSED_DIR
    / "ml_epigenomic_feature_table.csv"
)

BASELINE_FILE = (
    PROCESSED_DIR
    / "ml_feature_table.csv"
)

MODEL_FILE = (
    MODELS_DIR
    / "epigenomic_cox_model.pkl"
)

COEFFICIENT_FILE = (
    TABLES_DIR
    / "cox_model_coefficients.csv"
)

METRICS_FILE = (
    TABLES_DIR
    / "cox_model_metrics.csv"
)

PATIENT_RISK_FILE = (
    TABLES_DIR
    / "cox_patient_risk_scores.csv"
)

SELECTED_FEATURES_FILE = (
    TABLES_DIR
    / "cox_selected_features.csv"
)

KM_FIGURE = (
    FIGURES_DIR
    / "cox_risk_groups_survival.png"
)

TOP_FEATURES_FIGURE = (
    FIGURES_DIR
    / "cox_top_hazard_ratios.png"
)

REPORT_FILE = (
    REPORTS_DIR
    / "survival_model_report.md"
)


# =====================================
# SETTINGS
# =====================================

PATIENT_COLUMN = "patient_id"
EVENT_COLUMN = "outcome"

TEST_SIZE = 0.20
RANDOM_STATE = 42

TOP_METHYLATION_FEATURES = 50

COX_PENALIZER = 1.0

MAX_CATEGORICAL_LEVELS = 20

MIN_VARIANCE = 1e-8


# =====================================
# BASIC DATA FUNCTIONS
# =====================================

def load_dataset(
    file_path: Path,
) -> pd.DataFrame:
    """Load a CSV dataset."""

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
            f"Column '{PATIENT_COLUMN}' "
            "was not found."
        )

    df = df.copy()

    df[PATIENT_COLUMN] = (
        df[PATIENT_COLUMN]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return df


# =====================================
# OUTCOME
# =====================================

def add_outcome_column(
    epigenomic: pd.DataFrame,
    baseline: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add outcome from the Month 1 feature table
    to the epigenomic dataset.
    """

    if EVENT_COLUMN in epigenomic.columns:
        print(
            "Outcome column already exists."
        )

        return epigenomic

    print(
        "Outcome column not found."
    )

    print(
        "Adding outcome from Month 1 data..."
    )

    baseline = standardize_patient_ids(
        baseline
    )

    epigenomic = standardize_patient_ids(
        epigenomic
    )

    if EVENT_COLUMN not in baseline.columns:
        raise KeyError(
            f"Column '{EVENT_COLUMN}' "
            "was not found in ml_feature_table.csv."
        )

    outcome_table = baseline[
        [
            PATIENT_COLUMN,
            EVENT_COLUMN,
        ]
    ].copy()

    outcome_table[EVENT_COLUMN] = (
        pd.to_numeric(
            outcome_table[EVENT_COLUMN],
            errors="coerce",
        )
    )

    outcome_table = (
        outcome_table
        .dropna(
            subset=[
                PATIENT_COLUMN,
                EVENT_COLUMN,
            ]
        )
        .drop_duplicates(
            subset=[
                PATIENT_COLUMN
            ],
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
        epigenomic[
            EVENT_COLUMN
        ]
        .isna()
        .sum()
    )

    print(
        f"Patients without outcome: "
        f"{missing_outcomes:,}"
    )

    return epigenomic


# =====================================
# SURVIVAL TARGET
# =====================================

def find_survival_time_column(
    df: pd.DataFrame,
) -> str:
    """Find overall-survival duration column."""

    candidates = [
        "os_months",
        "overall_survival_months",
        "overall_survival",
        "survival_months",
    ]

    for column in candidates:
        if column in df.columns:
            return column

    raise KeyError(
        "No survival-time column found. "
        f"Checked: {candidates}"
    )


def validate_dataset(
    df: pd.DataFrame,
    duration_column: str,
) -> None:
    """Validate survival-analysis columns."""

    required = {
        PATIENT_COLUMN,
        EVENT_COLUMN,
        duration_column,
    }

    missing = required.difference(
        df.columns
    )

    if missing:
        raise KeyError(
            "Missing required columns: "
            f"{sorted(missing)}"
        )


def clean_survival_targets(
    df: pd.DataFrame,
    duration_column: str,
) -> pd.DataFrame:
    """Clean event and survival-time variables."""

    df = standardize_patient_ids(
        df
    )

    df[duration_column] = (
        pd.to_numeric(
            df[duration_column],
            errors="coerce",
        )
    )

    df[EVENT_COLUMN] = (
        pd.to_numeric(
            df[EVENT_COLUMN],
            errors="coerce",
        )
    )

    df = df.dropna(
        subset=[
            PATIENT_COLUMN,
            duration_column,
            EVENT_COLUMN,
        ]
    ).copy()

    df = df[
        df[duration_column] > 0
    ].copy()

    df[EVENT_COLUMN] = (
        df[EVENT_COLUMN]
        .astype(int)
    )

    df = df[
        df[EVENT_COLUMN].isin(
            [0, 1]
        )
    ].copy()

    df = (
        df.drop_duplicates(
            subset=[
                PATIENT_COLUMN
            ],
            keep="first",
        )
        .reset_index(drop=True)
    )

    return df


# =====================================
# METHYLATION FEATURES
# =====================================

def identify_methylation_columns(
    df: pd.DataFrame,
) -> list[str]:
    """Identify methylation columns."""

    return [
        column
        for column in df.columns
        if str(column).startswith(
            "meth_"
        )
    ]


def select_methylation_features(
    train_df: pd.DataFrame,
    methylation_columns: list[str],
    top_n: int,
) -> list[str]:
    """
    Select high-variance methylation features.

    Selection is performed using training data only.
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

    selected = (
        variance
        .sort_values(
            ascending=False
        )
        .head(top_n)
        .index
        .tolist()
    )

    return selected


# =====================================
# CLINICAL / GENOMIC FEATURE SELECTION
# =====================================

def identify_baseline_features(
    df: pd.DataFrame,
    duration_column: str,
) -> list[str]:
    """
    Select useful clinical and genomic predictors.

    Remove:
    - patient/sample identifiers
    - survival leakage variables
    - methylation variables
    """

    excluded_exact = {
        PATIENT_COLUMN,
        EVENT_COLUMN,
        duration_column,
        "sample_id",
        "case_id",
        "study_id",
        "project_id",
        "patient_barcode",
        "sample_barcode",
        "submitter_id",
        "tissue_source_site_code",
        "os_status",
        "vital_status",
        "survival_status",
        "os_status_patient",
    }

    excluded_fragments = [
        "sample_id",
        "case_id",
        "patient_id",
        "barcode",
        "submitter_id",
        "project_id",
        "study_id",
        "tissue_source_site",
    ]

    features = []

    for column in df.columns:

        name = str(column)

        lower_name = (
            name.lower()
        )

        if column in excluded_exact:
            continue

        if name.startswith(
            "meth_"
        ):
            continue

        if any(
            fragment in lower_name
            for fragment in excluded_fragments
        ):
            continue

        # Survival leakage
        if lower_name.endswith(
            "_months"
        ):
            continue

        if "survival" in lower_name:
            continue

        if lower_name.startswith(
            "os_"
        ):
            continue

        if lower_name.startswith(
            "dfs_"
        ):
            continue

        if lower_name.startswith(
            "pfs_"
        ):
            continue

        if lower_name.startswith(
            "dss_"
        ):
            continue

        features.append(
            column
        )

    return features


# =====================================
# PREPROCESSING
# =====================================

def prepare_train_test_features(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_columns: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """
    Prepare train/test Cox predictor matrices.

    All preprocessing decisions are based
    only on the training dataset.
    """

    train = train_df[
        feature_columns
    ].copy()

    test = test_df[
        feature_columns
    ].copy()

    numeric_columns = (
        train
        .select_dtypes(
            include=[
                np.number,
                "bool",
            ]
        )
        .columns
        .tolist()
    )

    categorical_columns = [
        column
        for column in train.columns
        if column not in numeric_columns
    ]

    # =================================
    # HIGH-CARDINALITY CATEGORICAL DATA
    # =================================

    high_cardinality = [
        column
        for column in categorical_columns
        if train[column].nunique(
            dropna=True
        ) > MAX_CATEGORICAL_LEVELS
    ]

    if high_cardinality:

        print(
            "\nDropping high-cardinality "
            "categorical columns:"
        )

        for column in high_cardinality:
            print(
                f"  - {column}"
            )

        train = train.drop(
            columns=high_cardinality
        )

        test = test.drop(
            columns=high_cardinality
        )

        categorical_columns = [
            column
            for column in categorical_columns
            if column not in high_cardinality
        ]

    # =================================
    # NUMERIC VARIABLES
    # =================================

    numeric_columns = [
        column
        for column in numeric_columns
        if column in train.columns
    ]

    for column in numeric_columns:

        train[column] = pd.to_numeric(
            train[column],
            errors="coerce",
        )

        test[column] = pd.to_numeric(
            test[column],
            errors="coerce",
        )

        median = train[
            column
        ].median()

        if pd.isna(median):
            median = 0.0

        train[column] = (
            train[column]
            .fillna(median)
        )

        test[column] = (
            test[column]
            .fillna(median)
        )

    # =================================
    # CATEGORICAL VARIABLES
    # =================================

    for column in categorical_columns:

        train[column] = (
            train[column]
            .astype("string")
        )

        test[column] = (
            test[column]
            .astype("string")
        )

        mode = train[
            column
        ].mode()

        if len(mode) > 0:
            fill_value = str(
                mode.iloc[0]
            )
        else:
            fill_value = "Missing"

        train[column] = (
            train[column]
            .fillna(fill_value)
        )

        test[column] = (
            test[column]
            .fillna(fill_value)
        )

    # =================================
    # ONE-HOT ENCODING
    # =================================

    if categorical_columns:

        train = pd.get_dummies(
            train,
            columns=categorical_columns,
            drop_first=True,
            dtype=float,
        )

        test = pd.get_dummies(
            test,
            columns=categorical_columns,
            drop_first=True,
            dtype=float,
        )

        test = test.reindex(
            columns=train.columns,
            fill_value=0.0,
        )

    # =================================
    # NUMERIC SAFETY
    # =================================

    train = train.apply(
        pd.to_numeric,
        errors="coerce",
    )

    test = test.apply(
        pd.to_numeric,
        errors="coerce",
    )

    train = train.replace(
        [
            np.inf,
            -np.inf,
        ],
        np.nan,
    )

    test = test.replace(
        [
            np.inf,
            -np.inf,
        ],
        np.nan,
    )

    train_medians = (
        train.median()
    )

    train = train.fillna(
        train_medians
    )

    test = test.fillna(
        train_medians
    )

    train = train.fillna(
        0.0
    )

    test = test.fillna(
        0.0
    )

    # =================================
    # LOW VARIANCE FILTER
    # =================================

    variance = train.var(
        axis=0
    )

    valid_columns = variance[
        variance > MIN_VARIANCE
    ].index.tolist()

    removed_count = (
        train.shape[1]
        - len(valid_columns)
    )

    if removed_count > 0:

        print(
            f"\nRemoved low-variance "
            f"predictors: {removed_count:,}"
        )

    train = train[
        valid_columns
    ].copy()

    test = test[
        valid_columns
    ].copy()

    print(
        f"\nPredictors after QC: "
        f"{train.shape[1]:,}"
    )

    return (
        train,
        test,
    )


# =====================================
# SURVIVAL DATASET
# =====================================

def build_survival_dataset(
    targets: pd.DataFrame,
    features: pd.DataFrame,
    duration_column: str,
) -> pd.DataFrame:
    """Combine predictors and survival outcomes."""

    survival_data = (
        features.copy()
    )

    survival_data[
        duration_column
    ] = (
        targets[
            duration_column
        ]
        .to_numpy()
    )

    survival_data[
        EVENT_COLUMN
    ] = (
        targets[
            EVENT_COLUMN
        ]
        .to_numpy()
    )

    return survival_data


# =====================================
# COX MODEL
# =====================================

def train_cox_model(
    train_data: pd.DataFrame,
    duration_column: str,
) -> CoxPHFitter:
    """Train penalized Cox proportional hazards model."""

    model = CoxPHFitter(
        penalizer=COX_PENALIZER
    )

    model.fit(
        train_data,
        duration_col=duration_column,
        event_col=EVENT_COLUMN,
        show_progress=False,
    )

    return model


def calculate_c_index(
    model: CoxPHFitter,
    data: pd.DataFrame,
    duration_column: str,
) -> float:
    """Calculate concordance index."""

    predictors = data.drop(
        columns=[
            duration_column,
            EVENT_COLUMN,
        ]
    )

    risk_scores = (
        model
        .predict_partial_hazard(
            predictors
        )
        .to_numpy()
        .ravel()
    )

    score = concordance_index(
        data[
            duration_column
        ],
        -risk_scores,
        data[
            EVENT_COLUMN
        ],
    )

    return float(
        score
    )


# =====================================
# COEFFICIENTS
# =====================================

def create_coefficient_table(
    model: CoxPHFitter,
) -> pd.DataFrame:
    """Create coefficient and hazard-ratio table."""

    summary = (
        model.summary
        .reset_index()
    )

    if "covariate" in summary.columns:

        summary = summary.rename(
            columns={
                "covariate":
                "feature"
            }
        )

    if "exp(coef)" in summary.columns:

        summary = summary.rename(
            columns={
                "exp(coef)":
                "hazard_ratio"
            }
        )

    if "feature" not in summary.columns:

        summary = summary.rename(
            columns={
                summary.columns[0]:
                "feature"
            }
        )

    keep_columns = [
        column
        for column in [
            "feature",
            "coef",
            "hazard_ratio",
            "se(coef)",
            "p",
        ]
        if column in summary.columns
    ]

    summary = summary[
        keep_columns
    ].copy()

    summary[
        "abs_coefficient"
    ] = (
        summary[
            "coef"
        ].abs()
    )

    summary = (
        summary
        .sort_values(
            "abs_coefficient",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )

    summary[
        "rank"
    ] = (
        summary.index
        + 1
    )

    return summary


# =====================================
# HAZARD-RATIO FIGURE
# =====================================

def save_hazard_ratio_plot(
    coefficients: pd.DataFrame,
) -> None:
    """Save hazard-ratio figure."""

    if (
        "hazard_ratio"
        not in coefficients.columns
    ):
        return

    plot_data = (
        coefficients
        .head(20)
        .sort_values(
            "hazard_ratio",
            ascending=True,
        )
    )

    plt.figure(
        figsize=(
            10,
            8,
        )
    )

    plt.barh(
        plot_data[
            "feature"
        ],
        plot_data[
            "hazard_ratio"
        ],
    )

    plt.axvline(
        1.0,
        linestyle="--",
    )

    plt.xlabel(
        "Hazard ratio"
    )

    plt.title(
        "Top Cox model hazard ratios"
    )

    plt.tight_layout()

    plt.savefig(
        TOP_FEATURES_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# =====================================
# RISK SCORES
# =====================================

def create_risk_scores(
    model: CoxPHFitter,
    test_data: pd.DataFrame,
    patient_ids: pd.Series,
    duration_column: str,
) -> pd.DataFrame:
    """Create patient-level Cox risk scores."""

    predictors = test_data.drop(
        columns=[
            duration_column,
            EVENT_COLUMN,
        ]
    )

    risk_scores = (
        model
        .predict_partial_hazard(
            predictors
        )
        .to_numpy()
        .ravel()
    )

    risk_table = pd.DataFrame(
        {
            PATIENT_COLUMN: (
                patient_ids
                .reset_index(
                    drop=True
                )
            ),
            duration_column: (
                test_data[
                    duration_column
                ]
                .reset_index(
                    drop=True
                )
            ),
            EVENT_COLUMN: (
                test_data[
                    EVENT_COLUMN
                ]
                .reset_index(
                    drop=True
                )
            ),
            "risk_score":
                risk_scores,
        }
    )

    median_risk = (
        risk_table[
            "risk_score"
        ]
        .median()
    )

    risk_table[
        "risk_group"
    ] = np.where(
        risk_table[
            "risk_score"
        ] >= median_risk,
        "High risk",
        "Low risk",
    )

    return risk_table


# =====================================
# KAPLAN-MEIER
# =====================================

def plot_kaplan_meier(
    risk_table: pd.DataFrame,
    duration_column: str,
) -> None:
    """Plot Kaplan-Meier curves by predicted risk."""

    plt.figure(
        figsize=(
            8,
            6,
        )
    )

    kmf = KaplanMeierFitter()

    for risk_group in [
        "Low risk",
        "High risk",
    ]:

        group = risk_table[
            risk_table[
                "risk_group"
            ] == risk_group
        ]

        if group.empty:
            continue

        kmf.fit(
            durations=group[
                duration_column
            ],
            event_observed=group[
                EVENT_COLUMN
            ],
            label=risk_group,
        )

        kmf.plot_survival_function(
            ci_show=True
        )

    plt.xlabel(
        "Overall survival time (months)"
    )

    plt.ylabel(
        "Survival probability"
    )

    plt.title(
        "Kaplan-Meier curves by predicted Cox risk"
    )

    plt.tight_layout()

    plt.savefig(
        KM_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# =====================================
# REPORT
# =====================================

def create_report(
    patient_count: int,
    train_count: int,
    test_count: int,
    event_count: int,
    censored_count: int,
    feature_count: int,
    methylation_feature_count: int,
    train_c_index: float,
    test_c_index: float,
    duration_column: str,
) -> str:
    """Create Markdown survival-model report."""

    report_lines = [
        "# Epigenomic Survival Model Report",
        "",
        "## Cohort",
        "",
        (
            f"- Total patients: "
            f"{patient_count:,}"
        ),
        (
            f"- Training patients: "
            f"{train_count:,}"
        ),
        (
            f"- Testing patients: "
            f"{test_count:,}"
        ),
        (
            f"- Events/deaths: "
            f"{event_count:,}"
        ),
        (
            f"- Censored/living: "
            f"{censored_count:,}"
        ),
        "",
        "## Survival outcome",
        "",
        (
            "- Duration variable: "
            f"`{duration_column}`"
        ),
        (
            "- Event variable: "
            f"`{EVENT_COLUMN}`"
        ),
        "- Event = 1: deceased",
        "- Event = 0: censored/living",
        "",
        "## Predictors",
        "",
        (
            "- Final encoded predictors: "
            f"{feature_count:,}"
        ),
        (
            "- Selected methylation features: "
            f"{methylation_feature_count:,}"
        ),
        "",
        "## Model",
        "",
        (
            "- Algorithm: Cox proportional "
            "hazards regression"
        ),
        (
            f"- Penalizer: "
            f"{COX_PENALIZER}"
        ),
        "",
        "## Performance",
        "",
        (
            "- Training C-index: "
            f"{train_c_index:.3f}"
        ),
        (
            "- Testing C-index: "
            f"{test_c_index:.3f}"
        ),
        "",
        "## Interpretation",
        "",
        (
            "A concordance index close to 0.5 "
            "indicates approximately random risk "
            "ranking. Higher values indicate better "
            "discrimination between patients with "
            "shorter and longer survival."
        ),
        "",
        (
            "Hazard ratios greater than 1 indicate "
            "increased estimated hazard, whereas "
            "hazard ratios below 1 indicate decreased "
            "estimated hazard."
        ),
        "",
        "## Leakage controls",
        "",
        (
            "- Patient and sample identifiers were "
            "excluded from model predictors."
        ),
        (
            "- Other survival-related variables were "
            "excluded from predictors."
        ),
        (
            "- High-cardinality categorical variables "
            "were removed."
        ),
        (
            "- Methylation feature selection was "
            "performed using training data only."
        ),
        (
            "- Missing-value imputation was based on "
            "training data only."
        ),
        "",
    ]

    return "\n".join(
        report_lines
    )


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Train epigenomic Cox survival model."""

    print(
        "Loading epigenomic feature table..."
    )

    dataset = load_dataset(
        INPUT_FILE
    )

    print(
        "Loading Month 1 outcome data..."
    )

    baseline = load_dataset(
        BASELINE_FILE
    )

    dataset = add_outcome_column(
        epigenomic=dataset,
        baseline=baseline,
    )

    duration_column = (
        find_survival_time_column(
            dataset
        )
    )

    print(
        f"Survival time column: "
        f"{duration_column}"
    )

    validate_dataset(
        dataset,
        duration_column,
    )

    dataset = clean_survival_targets(
        dataset,
        duration_column,
    )

    print(
        "\nSurvival cohort"
    )

    print(
        "-" * 40
    )

    print(
        f"Patients with usable survival data: "
        f"{len(dataset):,}"
    )

    event_count = int(
        dataset[
            EVENT_COLUMN
        ]
        .sum()
    )

    censored_count = int(
        len(dataset)
        - event_count
    )

    print(
        f"Events/deaths: "
        f"{event_count:,}"
    )

    print(
        f"Censored/living: "
        f"{censored_count:,}"
    )

    # =================================
    # TRAIN / TEST SPLIT FIRST
    # =================================

    indices = np.arange(
        len(dataset)
    )

    (
        train_indices,
        test_indices,
    ) = train_test_split(
        indices,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=dataset[
            EVENT_COLUMN
        ],
    )

    train_df = (
        dataset
        .iloc[
            train_indices
        ]
        .reset_index(
            drop=True
        )
    )

    test_df = (
        dataset
        .iloc[
            test_indices
        ]
        .reset_index(
            drop=True
        )
    )

    print(
        f"\nTraining patients: "
        f"{len(train_df):,}"
    )

    print(
        f"Testing patients: "
        f"{len(test_df):,}"
    )

    # =================================
    # FEATURE SELECTION
    # =================================

    methylation_columns = (
        identify_methylation_columns(
            dataset
        )
    )

    print(
        f"\nAvailable methylation features: "
        f"{len(methylation_columns):,}"
    )

    selected_methylation = (
        select_methylation_features(
            train_df=train_df,
            methylation_columns=(
                methylation_columns
            ),
            top_n=(
                TOP_METHYLATION_FEATURES
            ),
        )
    )

    print(
        f"Selected methylation features: "
        f"{len(selected_methylation):,}"
    )

    baseline_features = (
        identify_baseline_features(
            dataset,
            duration_column,
        )
    )

    print(
        f"Candidate baseline features: "
        f"{len(baseline_features):,}"
    )

    feature_columns = (
        baseline_features
        + selected_methylation
    )

    feature_columns = list(
        dict.fromkeys(
            feature_columns
        )
    )

    # =================================
    # PREPROCESSING
    # =================================

    (
        train_features,
        test_features,
    ) = prepare_train_test_features(
        train_df=train_df,
        test_df=test_df,
        feature_columns=feature_columns,
    )

    selected_features_table = (
        pd.DataFrame(
            {
                "feature":
                    train_features.columns
            }
        )
    )

    selected_features_table[
        "is_methylation"
    ] = (
        selected_features_table[
            "feature"
        ]
        .str.startswith(
            "meth_"
        )
    )

    selected_features_table.to_csv(
        SELECTED_FEATURES_FILE,
        index=False,
    )

    # =================================
    # SURVIVAL MATRICES
    # =================================

    train_data = build_survival_dataset(
        targets=train_df,
        features=train_features,
        duration_column=duration_column,
    )

    test_data = build_survival_dataset(
        targets=test_df,
        features=test_features,
        duration_column=duration_column,
    )

    # =================================
    # MODEL
    # =================================

    print(
        "\nTraining Cox proportional "
        "hazards model..."
    )

    model = train_cox_model(
        train_data=train_data,
        duration_column=duration_column,
    )

    # =================================
    # PERFORMANCE
    # =================================

    train_c_index = calculate_c_index(
        model=model,
        data=train_data,
        duration_column=duration_column,
    )

    test_c_index = calculate_c_index(
        model=model,
        data=test_data,
        duration_column=duration_column,
    )

    print(
        f"\nTraining C-index: "
        f"{train_c_index:.3f}"
    )

    print(
        f"Testing C-index: "
        f"{test_c_index:.3f}"
    )

    # =================================
    # COEFFICIENTS
    # =================================

    coefficients = (
        create_coefficient_table(
            model
        )
    )

    coefficients.to_csv(
        COEFFICIENT_FILE,
        index=False,
    )

    save_hazard_ratio_plot(
        coefficients
    )

    # =================================
    # METRICS
    # =================================

    metrics = pd.DataFrame(
        {
            "metric": [
                "training_c_index",
                "testing_c_index",
                "patients",
                "training_patients",
                "testing_patients",
                "events",
                "censored",
                "final_predictors",
                "selected_methylation_features",
            ],
            "value": [
                train_c_index,
                test_c_index,
                len(dataset),
                len(train_data),
                len(test_data),
                event_count,
                censored_count,
                train_features.shape[1],
                len(
                    selected_methylation
                ),
            ],
        }
    )

    metrics.to_csv(
        METRICS_FILE,
        index=False,
    )

    # =================================
    # RISK SCORES
    # =================================

    risk_table = create_risk_scores(
        model=model,
        test_data=test_data,
        patient_ids=test_df[
            PATIENT_COLUMN
        ],
        duration_column=duration_column,
    )

    risk_table.to_csv(
        PATIENT_RISK_FILE,
        index=False,
    )

    # =================================
    # KAPLAN-MEIER
    # =================================

    plot_kaplan_meier(
        risk_table=risk_table,
        duration_column=duration_column,
    )

    # =================================
    # SAVE MODEL
    # =================================

    model_package = {
        "model": model,
        "duration_column":
            duration_column,
        "event_column":
            EVENT_COLUMN,
        "predictor_columns":
            train_features.columns.tolist(),
        "selected_methylation":
            selected_methylation,
        "cox_penalizer":
            COX_PENALIZER,
    }

    joblib.dump(
        model_package,
        MODEL_FILE,
    )

    # =================================
    # REPORT
    # =================================

    report = create_report(
        patient_count=len(dataset),
        train_count=len(train_data),
        test_count=len(test_data),
        event_count=event_count,
        censored_count=censored_count,
        feature_count=(
            train_features.shape[1]
        ),
        methylation_feature_count=(
            len(selected_methylation)
        ),
        train_c_index=train_c_index,
        test_c_index=test_c_index,
        duration_column=duration_column,
    )

    REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    # =================================
    # DONE
    # =================================

    print(
        "\nFinished successfully."
    )

    print(
        "\nCox model saved to:\n"
        f"{MODEL_FILE}"
    )

    print(
        "\nCoefficient table saved to:\n"
        f"{COEFFICIENT_FILE}"
    )

    print(
        "\nMetrics saved to:\n"
        f"{METRICS_FILE}"
    )

    print(
        "\nSelected features saved to:\n"
        f"{SELECTED_FEATURES_FILE}"
    )

    print(
        "\nRisk scores saved to:\n"
        f"{PATIENT_RISK_FILE}"
    )

    print(
        "\nKaplan-Meier figure saved to:\n"
        f"{KM_FIGURE}"
    )

    print(
        "\nHazard-ratio figure saved to:\n"
        f"{TOP_FEATURES_FIGURE}"
    )

    print(
        "\nTechnical report saved to:\n"
        f"{REPORT_FILE}"
    )


if __name__ == "__main__":
    main()