# Final Epigenomic Classification Model

## Model specification

- Algorithm: Random Forest classifier
- Probability calibration: isotonic regression
- Random Forest estimators: 300
- Calibration folds: 4

## Cohort

- Development patients: 707
- Positive outcomes: 88
- Negative outcomes: 619

## Predictors

- Final raw predictors: 546
- Selected methylation predictors: 500

## Development threshold

- Threshold was selected from leakage-aware out-of-fold calibrated probabilities.
- Selected probability threshold: 0.115

## Out-of-fold performance

- ROC AUC: 0.853
- PR AUC: 0.660
- Brier score: 0.0674
- Log loss: 0.2621
- Accuracy: 0.784
- Balanced accuracy: 0.779
- Precision: 0.338
- Recall: 0.773
- F1-score: 0.471

## Methodological safeguards

- Patient/sample identifiers were excluded from model predictors.
- Survival and direct outcome-proxy variables were excluded.
- Methylation selection for OOF predictions was performed within each training fold.
- Isotonic calibration was fitted without access to outer-test outcomes.

## Interpretation

This model is the final development-stage clinical-genomic-epigenomic classifier.

The probability threshold is a development-stage operating threshold derived from internal cross-validation. It should not be interpreted as an externally validated clinical decision threshold.

Independent external validation is required before any clinical interpretation or deployment.
