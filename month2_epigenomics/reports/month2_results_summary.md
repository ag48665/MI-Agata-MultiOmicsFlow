# Month 2 Results Summary

## Results

Integration of DNA methylation substantially improved predictive discrimination compared with the clinical-genomic baseline.

In leakage-aware five-fold cross-validation, the clinical-genomic model achieved a mean ROC AUC of 0.548 ± 0.026, whereas the clinical-genomic-epigenomic model achieved 0.866 ± 0.051.

This corresponded to an absolute ROC AUC increase of 0.318 after inclusion of methylation features.

## Classification threshold

At the conventional probability threshold of 0.50, the epigenomic classifier showed high overall accuracy but low sensitivity for the positive class.

Mean recall was 0.169, with balanced accuracy of 0.584.

Training-derived threshold optimisation improved minority-class detection.

The mean optimized uncalibrated threshold was 0.278 ± 0.024.

At the optimized threshold, recall increased to 0.774, balanced accuracy increased to 0.785, and F1-score increased to 0.509.

## Probability calibration

Probability calibration improved the quality of predicted risk estimates.

The uncalibrated model had a mean Brier score of 0.0922, whereas isotonic calibration reduced the mean Brier score to 0.0683.

Isotonic calibration achieved a mean ROC AUC of 0.855 and mean PR AUC of 0.671.

## Final development model

The final development-stage classifier consisted of a Random Forest model with isotonic probability calibration.

The final calibrated operating threshold was 0.115.

Leakage-aware out-of-fold performance of the final model included ROC AUC 0.853, PR AUC 0.660, Brier score 0.0674, balanced accuracy 0.779, precision 0.338, recall 0.773, and F1-score 0.471.

## Key results table

| result                                |   value |      std |
|:--------------------------------------|--------:|---------:|
| Clinical-genomic ROC AUC              |  0.5482 |   0.0264 |
| Clinical-genomic-epigenomic ROC AUC   |  0.8657 |   0.0505 |
| ROC AUC improvement after methylation |  0.3175 | nan      |
| Default threshold recall              |  0.1686 |   0.0763 |
| Optimized threshold recall            |  0.7739 |   0.0752 |
| Default threshold balanced accuracy   |  0.5835 |   0.0396 |
| Optimized threshold balanced accuracy |  0.7853 |   0.0655 |
| Optimized uncalibrated threshold      |  0.278  |   0.0239 |
| Uncalibrated Brier score              |  0.0922 |   0.0045 |
| Sigmoid Brier score                   |  0.0718 |   0.0023 |
| Isotonic Brier score                  |  0.0683 |   0.0036 |
| Final calibrated ROC AUC              |  0.853  | nan      |
| Final calibrated PR AUC               |  0.66   | nan      |
| Final calibrated Brier score          |  0.0674 | nan      |
| Final calibrated threshold            |  0.115  | nan      |
| Final balanced accuracy               |  0.7789 | nan      |
| Final precision                       |  0.3383 | nan      |
| Final recall                          |  0.7727 | nan      |
| Final F1-score                        |  0.4706 | nan      |

## Recommended figures

1. `model_validation_roc_auc.png`
2. `threshold_validation_roc_curve.png`
3. `threshold_validation_precision_recall.png`
4. `confusion_matrix_optimized_threshold.png`
5. `model_calibration_curve.png`
6. `cox_risk_groups_survival.png`
7. `cox_top_hazard_ratios.png`
8. `combined_feature_importance.png`

## Interpretation

The results support the hypothesis that DNA methylation contains predictive information beyond the clinical-genomic baseline.

The decline from the original single-split ROC AUC of 1.0 to approximately 0.86 under leakage-aware cross-validation demonstrates the importance of rigorous validation.

Threshold optimisation substantially improved positive-class sensitivity, while probability calibration improved agreement between predicted probabilities and observed outcomes.

## Limitations

All reported model-performance estimates are based on internal validation within the available cohort.

External validation on an independent dataset is required before assessing transportability or clinical utility.

The optimized classification threshold is a development-stage operating point and should not be interpreted as a validated clinical cutoff.
