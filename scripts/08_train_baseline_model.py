from pathlib import Path
import joblib
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


# ======================================
# PATHS
# ======================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "ml_feature_table.csv"
)

MODEL_DIR = PROJECT_ROOT / "outputs" / "models"
FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"
TABLE_DIR = PROJECT_ROOT / "outputs" / "tables"

MODEL_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


# ======================================
# LOAD DATA
# ======================================

print("Loading ML feature table...")

df = pd.read_csv(INPUT_FILE)

print(df.shape)


# ======================================
# PREPARE DATA
# ======================================

drop_columns = ["patient_id", "outcome"]

X = df.drop(columns=drop_columns)
y = df["outcome"]

from pandas.api.types import is_numeric_dtype


categorical_features = [
    column
    for column in X.columns
    if not is_numeric_dtype(X[column])
]

numeric_features = [
    column
    for column in X.columns
    if is_numeric_dtype(X[column])
]

print("\nCategorical:")
print(categorical_features)

print("\nNumeric:")
print(numeric_features)


numeric_transformer = Pipeline(
    steps=[
        ("imputer", SimpleImputer(strategy="median"))
    ]
)

categorical_transformer = Pipeline(
    steps=[
        ("imputer", SimpleImputer(strategy="most_frequent")),
        (
            "encoder",
            OneHotEncoder(
                handle_unknown="ignore"
            ),
        ),
    ]
)

preprocessor = ColumnTransformer(
    transformers=[
        (
            "num",
            numeric_transformer,
            numeric_features,
        ),
        (
            "cat",
            categorical_transformer,
            categorical_features,
        ),
    ]
)


X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y,
)

print("\nTraining samples:", len(X_train))
print("Testing samples:", len(X_test))


# ======================================
# MODEL
# ======================================

model = RandomForestClassifier(
    n_estimators=300,
    max_depth=None,
    random_state=42,
    class_weight="balanced",
)

pipeline = Pipeline(
    steps=[
        ("preprocessor", preprocessor),
        ("model", model),
    ]
)

print("\nTraining Random Forest...")

pipeline.fit(X_train, y_train)

print("Training complete.")

# ======================================
# PREDICTIONS
# ======================================

print("\nEvaluating model...")

y_pred = pipeline.predict(X_test)
y_prob = pipeline.predict_proba(X_test)[:, 1]

accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred, zero_division=0)
recall = recall_score(y_test, y_pred, zero_division=0)
f1 = f1_score(y_test, y_pred, zero_division=0)
roc = roc_auc_score(y_test, y_prob)

print("\nModel performance")
print(f"Accuracy : {accuracy:.3f}")
print(f"Precision: {precision:.3f}")
print(f"Recall   : {recall:.3f}")
print(f"F1-score : {f1:.3f}")
print(f"ROC AUC  : {roc:.3f}")


# ======================================
# SAVE METRICS
# ======================================

metrics = pd.DataFrame(
    {
        "metric": [
            "accuracy",
            "precision",
            "recall",
            "f1_score",
            "roc_auc",
        ],
        "value": [
            accuracy,
            precision,
            recall,
            f1,
            roc,
        ],
    }
)

metrics.to_csv(
    TABLE_DIR / "baseline_model_metrics.csv",
    index=False,
)


# ======================================
# CLASSIFICATION REPORT
# ======================================

report = classification_report(
    y_test,
    y_pred,
    output_dict=True,
    zero_division=0,
)

pd.DataFrame(report).transpose().to_csv(
    TABLE_DIR / "classification_report.csv"
)


# ======================================
# CONFUSION MATRIX
# ======================================

cm = confusion_matrix(y_test, y_pred)

disp = ConfusionMatrixDisplay(cm)

disp.plot()

plt.title("Random Forest Confusion Matrix")

plt.tight_layout()

plt.savefig(
    FIGURE_DIR / "confusion_matrix.png",
    dpi=300,
)

plt.close()


# ======================================
# FEATURE IMPORTANCE
# ======================================

feature_names = (
    pipeline.named_steps["preprocessor"]
    .get_feature_names_out()
)

importances = (
    pipeline.named_steps["model"]
    .feature_importances_
)

importance_df = pd.DataFrame(
    {
        "feature": feature_names,
        "importance": importances,
    }
).sort_values(
    "importance",
    ascending=False,
)

importance_df.to_csv(
    TABLE_DIR / "feature_importance.csv",
    index=False,
)

top20 = importance_df.head(20)

plt.figure(figsize=(10, 8))

plt.barh(
    top20["feature"][::-1],
    top20["importance"][::-1],
)

plt.xlabel("Importance")

plt.title("Top 20 Feature Importances")

plt.tight_layout()

plt.savefig(
    FIGURE_DIR / "feature_importance.png",
    dpi=300,
)

plt.close()


# ======================================
# SAVE MODEL
# ======================================

joblib.dump(
    pipeline,
    MODEL_DIR / "baseline_random_forest.pkl",
)


# ======================================
# FINISH
# ======================================

print("\nModel saved:")
print(MODEL_DIR / "baseline_random_forest.pkl")

print("\nOutputs created:")

print(TABLE_DIR / "baseline_model_metrics.csv")
print(TABLE_DIR / "classification_report.csv")
print(TABLE_DIR / "feature_importance.csv")

print(FIGURE_DIR / "confusion_matrix.png")
print(FIGURE_DIR / "feature_importance.png")

print("\nDone!")