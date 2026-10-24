# Results

## Predictive contribution of DNA methylation

Leakage-aware five-fold cross-validation demonstrated a substantial increase in discriminatory performance after incorporation of DNA methylation features. The clinical-genomic baseline model achieved a mean ROC AUC of 0.548 ± 0.026. In comparison, the clinical-genomic-epigenomic model achieved a mean ROC AUC of 0.866 ± 0.051. The absolute improvement in mean ROC AUC associated with inclusion of methylation features was therefore 0.318.

The clinical-genomic-epigenomic model also achieved a mean balanced accuracy of 0.584 ± 0.040. These results indicated that the methylation layer contributed predictive information beyond that available from the clinical-genomic baseline.

## Effect of classification-threshold optimisation

At the default classification threshold of 0.50, the epigenomic model achieved a mean accuracy of 0.895, but minority-class sensitivity remained low. Mean recall at the default threshold was 0.169 ± 0.076, while mean balanced accuracy was 0.584 ± 0.040 and mean F1-score was 0.281 ± 0.119.

Threshold optimisation performed within the training folds identified a mean operating threshold of 0.278 ± 0.024. At this threshold, mean recall increased to 0.774 ± 0.075, mean balanced accuracy increased to 0.785 ± 0.066, and mean F1-score increased to 0.509 ± 0.132. Precision decreased from 0.900 at the default threshold to 0.400 after optimisation, demonstrating the expected trade-off between sensitivity and positive predictive value.

ROC AUC remained 0.866, and PR AUC remained 0.669, as these measures are not determined by a single classification threshold.

## Probability calibration

Probability calibration improved probabilistic model performance. The uncalibrated Random Forest produced a mean Brier score of 0.0922 ± 0.0045 and a mean log loss of 0.3347 ± 0.0115.

Isotonic calibration produced the lowest mean Brier score, decreasing it to 0.0683 ± 0.0036. Mean log loss was also reduced to 0.2510 ± 0.0215. Following isotonic calibration, mean ROC AUC was 0.855 ± 0.053 and mean PR AUC was 0.671 ± 0.040.

On the basis of the Brier score, isotonic regression was selected as the final probability-calibration method.

## Final development-stage model

The final model consisted of a Random Forest classifier combined with isotonic probability calibration. A calibrated operating threshold of 0.115 was selected using leakage-aware out-of-fold predictions.

The final out-of-fold ROC AUC was 0.853, with a PR AUC of 0.660. The final Brier score was 0.0674 and log loss was 0.2621.

At the selected calibrated threshold, the model achieved an accuracy of 0.784, balanced accuracy of 0.779, precision of 0.338, recall of 0.773, and F1-score of 0.471.

Overall, the validation analyses showed that the initially observed near-perfect single-split performance was reduced under more stringent leakage-aware cross-validation, but the epigenomic model retained substantially stronger discrimination than the clinical-genomic baseline.