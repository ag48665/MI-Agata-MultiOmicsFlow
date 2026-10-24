# Classification Threshold Optimisation Report

## Objective

The leakage-aware validation in Step 21 showed that the epigenomic model achieved strong ROC AUC but relatively low recall when the default classification threshold of 0.50 was used.

This analysis evaluates whether a training-derived decision threshold improves classification of the minority outcome class.

## Validation design

- Patients: 707
- Outer validation folds: 5
- Inner threshold-tuning folds: 4
- Methylation feature selection was performed using outer-training data only.
- Threshold optimisation was performed using predictions generated only from training data.
- The optimized threshold was then evaluated on the independent outer test fold.

## Outcome distribution

- Positive-class prevalence: 0.124

## Cross-validation results

| strategy     |   folds |   accuracy_mean |   accuracy_std |   balanced_accuracy_mean |   balanced_accuracy_std |   precision_mean |   precision_std |   recall_mean |   recall_std |   f1_score_mean |   f1_score_std |   roc_auc_mean |   roc_auc_std |   pr_auc_mean |   pr_auc_std |   threshold_mean |   threshold_std |
|:-------------|--------:|----------------:|---------------:|-------------------------:|------------------------:|-----------------:|----------------:|--------------:|-------------:|----------------:|---------------:|---------------:|--------------:|--------------:|-------------:|-----------------:|----------------:|
| default_0.50 |       5 |          0.8953 |         0.0095 |                   0.5835 |                  0.0396 |           0.9    |          0.2236 |        0.1686 |       0.0763 |          0.2813 |         0.1187 |         0.8657 |        0.0505 |         0.669 |        0.071 |            0.5   |          0      |
| optimized    |       5 |          0.7935 |         0.1044 |                   0.7853 |                  0.0655 |           0.4003 |          0.1613 |        0.7739 |       0.0752 |          0.5094 |         0.1322 |         0.8657 |        0.0505 |         0.669 |        0.071 |            0.278 |          0.0239 |

## Fold-level results

|   fold | strategy     |   threshold |   accuracy |   balanced_accuracy |   precision |   recall |   f1_score |   roc_auc |   pr_auc |
|-------:|:-------------|------------:|-----------:|--------------------:|------------:|---------:|-----------:|----------:|---------:|
|      1 | default_0.50 |        0.5  |     0.9014 |              0.6111 |      1      |   0.2222 |     0.3636 |    0.9444 |   0.7614 |
|      1 | optimized    |        0.29 |     0.8873 |              0.8642 |      0.5357 |   0.8333 |     0.6522 |    0.9444 |   0.7614 |
|      2 | default_0.50 |        0.5  |     0.9014 |              0.6111 |      1      |   0.2222 |     0.3636 |    0.8058 |   0.6347 |
|      2 | optimized    |        0.25 |     0.662  |              0.6877 |      0.2321 |   0.7222 |     0.3514 |    0.8058 |   0.6347 |
|      3 | default_0.50 |        0.5  |     0.8936 |              0.5588 |      1      |   0.1176 |     0.2105 |    0.8712 |   0.6542 |
|      3 | optimized    |        0.28 |     0.8014 |              0.811  |      0.359  |   0.8235 |     0.5    |    0.8712 |   0.6542 |
|      4 | default_0.50 |        0.5  |     0.8794 |              0.5254 |      0.5    |   0.0588 |     0.1053 |    0.8468 |   0.5793 |
|      4 | optimized    |        0.26 |     0.7163 |              0.7626 |      0.2745 |   0.8235 |     0.4118 |    0.8468 |   0.5793 |
|      5 | default_0.50 |        0.5  |     0.9007 |              0.6111 |      1      |   0.2222 |     0.3636 |    0.8604 |   0.7157 |
|      5 | optimized    |        0.31 |     0.9007 |              0.8008 |      0.6    |   0.6667 |     0.6316 |    0.8604 |   0.7157 |

## Probability discrimination

- Pooled outer-fold ROC AUC: 0.865
- Pooled outer-fold PR AUC: 0.667

## Optimized threshold

- Mean training-derived optimized threshold: 0.278

## Interpretation

ROC AUC describes ranking performance and does not depend on a single classification threshold.

Because the outcome classes are imbalanced, accuracy alone can be misleading. Balanced accuracy, recall, F1-score and precision-recall performance should therefore be considered alongside ROC AUC.

The optimized threshold should only be preferred if it improves clinically relevant metrics such as recall and balanced accuracy on the independent outer folds.

The threshold estimated here should be considered a validation-stage threshold rather than a final clinical decision threshold.
