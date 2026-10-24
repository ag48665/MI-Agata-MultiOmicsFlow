# Model Probability Calibration Report

## Objective

This analysis evaluates whether probability calibration improves the reliability of predicted outcome probabilities produced by the clinical-genomic-epigenomic Random Forest model.

## Cohort

- Patients: 707
- Positive outcomes: 88
- Negative outcomes: 619

## Validation design

- Outer validation folds: 5
- Calibration folds: 4
- Methylation feature selection was performed using outer-training data only.
- Calibration was fitted without access to outer-test outcomes.

## Calibration methods

- Uncalibrated Random Forest
- Sigmoid / Platt-style calibration
- Isotonic regression calibration

## Cross-validation summary

| method       |   folds |   roc_auc_mean |   roc_auc_std |   pr_auc_mean |   pr_auc_std |   brier_score_mean |   brier_score_std |   log_loss_mean |   log_loss_std |
|:-------------|--------:|---------------:|--------------:|--------------:|-------------:|-------------------:|------------------:|----------------:|---------------:|
| uncalibrated |       5 |         0.8657 |        0.0505 |        0.669  |       0.071  |             0.0922 |            0.0045 |          0.3347 |         0.0115 |
| sigmoid      |       5 |         0.8632 |        0.0478 |        0.6694 |       0.0395 |             0.0718 |            0.0023 |          0.2587 |         0.0137 |
| isotonic     |       5 |         0.8547 |        0.053  |        0.6709 |       0.0405 |             0.0683 |            0.0036 |          0.251  |         0.0215 |

## Fold-level results

|   roc_auc |   pr_auc |   brier_score |   log_loss |   fold | method       |
|----------:|---------:|--------------:|-----------:|-------:|:-------------|
|    0.9444 |   0.7614 |        0.088  |     0.3245 |      1 | uncalibrated |
|    0.9167 |   0.7213 |        0.0704 |     0.2501 |      1 | sigmoid      |
|    0.8958 |   0.7218 |        0.0655 |     0.2367 |      1 | isotonic     |
|    0.8058 |   0.6347 |        0.0974 |     0.3482 |      2 | uncalibrated |
|    0.7854 |   0.6455 |        0.0756 |     0.2831 |      2 | sigmoid      |
|    0.7634 |   0.622  |        0.074  |     0.2853 |      2 | isotonic     |
|    0.8712 |   0.6542 |        0.0915 |     0.3312 |      3 | uncalibrated |
|    0.87   |   0.6383 |        0.0722 |     0.2538 |      3 | sigmoid      |
|    0.8871 |   0.6822 |        0.0655 |     0.2315 |      3 | isotonic     |
|    0.8468 |   0.5793 |        0.0962 |     0.3456 |      4 | uncalibrated |
|    0.8686 |   0.6393 |        0.0712 |     0.2543 |      4 | sigmoid      |
|    0.8665 |   0.6382 |        0.0698 |     0.2571 |      4 | isotonic     |
|    0.8604 |   0.7157 |        0.0877 |     0.324  |      5 | uncalibrated |
|    0.8753 |   0.7027 |        0.0698 |     0.2522 |      5 | sigmoid      |
|    0.8604 |   0.6902 |        0.0668 |     0.2442 |      5 | isotonic     |

## Best calibration results

- Lowest mean Brier score: `isotonic` (0.0683)
- Lowest mean log loss: `isotonic` (0.2510)

## Interpretation

Lower Brier score and lower log loss indicate better probabilistic predictions.

ROC AUC and PR AUC primarily describe discrimination, whereas calibration assesses whether predicted probabilities correspond to observed outcome frequencies.

A calibrated probability should not automatically be interpreted as an externally validated clinical risk estimate. Independent validation would still be required before clinical use.
