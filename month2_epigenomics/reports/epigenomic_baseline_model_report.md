# Epigenomic Baseline Model Report

## Cohort

- Matched patients: 707
- Training patients: 565
- Testing patients: 142
- Clinical-genomic input features: 17
- Epigenomic input features: 559

## Model settings

- Algorithm: Random Forest Classifier
- Trees: 500
- Class weight: balanced
- Train/test split: 80/20
- Shared patient split for both models

## Clinical-genomic model

- Accuracy: 0.810
- Balanced accuracy: 0.487
- Precision: 0.091
- Recall: 0.056
- F1-score: 0.069
- ROC AUC: 0.450

## Clinical-genomic-epigenomic model

- Accuracy: 0.986
- Balanced accuracy: 0.944
- Precision: 1.000
- Recall: 0.889
- F1-score: 0.941
- ROC AUC: 1.000

## Epigenomic contribution

- ROC AUC change after adding methylation: +0.550

## Interpretation note

The target is overall survival status. This is a baseline classification experiment and not a full time-to-event survival analysis.
