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

# Discussion

The integration of DNA methylation with clinical and genomic data resulted in a marked improvement in predictive discrimination. Under leakage-aware five-fold cross-validation, mean ROC AUC increased from 0.548 for the clinical-genomic model to 0.866 for the clinical-genomic-epigenomic model, corresponding to an absolute increase of 0.318. This finding supports the interpretation that epigenomic variation contains information that is complementary to the baseline clinical and genomic feature set.

An important methodological observation was the difference between the original single-split model performance and the results obtained under stricter cross-validation. The earlier ROC AUC of approximately 1.0 was not reproduced when feature selection, preprocessing and model evaluation were constrained within cross-validation folds. The reduction to a mean ROC AUC of approximately 0.866 illustrates the importance of leakage-aware validation in high-dimensional multi-omics analyses. Although the revised performance estimate was lower, it remained substantially above that of the clinical-genomic baseline.

Class imbalance also had a major effect on model interpretation. At the default probability threshold of 0.50, overall accuracy was high (0.895), but mean recall was only 0.169. This demonstrated that conventional accuracy alone would have overstated model usefulness for detecting the minority outcome class. Optimising the operating threshold increased recall to 0.774 and balanced accuracy to 0.785. The associated reduction in precision to 0.400 reflects the expected sensitivity-specificity trade-off and highlights the need to define the intended use case before selecting an operating threshold.

Probability calibration provided an additional improvement in model interpretability. Isotonic calibration reduced the mean Brier score from 0.0922 to 0.0683 and reduced mean log loss from 0.3347 to 0.2510. The small reduction in ROC AUC following isotonic calibration is not inconsistent with successful calibration, because calibration primarily concerns agreement between predicted probabilities and observed event frequencies rather than ranking performance.

The final development-stage model achieved an out-of-fold ROC AUC of 0.853 and PR AUC of 0.660, with a Brier score of 0.0674. At the final calibrated threshold of 0.115, balanced accuracy was 0.779 and recall was 0.773. These values suggest that the model retained meaningful discrimination while providing a more balanced operating point than the original default-threshold classifier.

Several limitations remain. First, all estimates were obtained through internal validation within the same overall cohort. Internal cross-validation reduces optimism but does not establish transportability to an independent population. Second, the high dimensionality of methylation data relative to the number of patients increases the risk of overfitting and instability in individual selected features. Third, the operating threshold was selected to maximise balanced accuracy within the development data and should not be interpreted as a clinically validated cutoff. Fourth, calibration performance may change when disease prevalence or patient characteristics differ in external cohorts.

The results therefore support further investigation of methylation-enhanced multi-omics prediction, but they do not establish clinical utility. Independent external validation should be the next major analytical step. In addition, future work could examine stability of selected methylation loci, biological pathway enrichment, alternative feature-selection strategies and model architectures, and validation across independent breast cancer cohorts.

# Conclusions

Integration of DNA methylation into the clinical-genomic workflow substantially improved predictive discrimination in the available development cohort. Leakage-aware cross-validation increased mean ROC AUC from 0.548 for the clinical-genomic baseline to 0.866 after inclusion of epigenomic features.

Threshold optimisation improved detection of the minority outcome class, while isotonic calibration improved the quality of predicted probabilities. The final calibrated model achieved an internal out-of-fold ROC AUC of 0.853, PR AUC of 0.660, balanced accuracy of 0.779, and recall of 0.773.

These findings support the potential value of DNA methylation as an additional predictive layer in multi-omics modelling. However, the model remains a development-stage research model. Independent external validation is required before the reported probabilities, selected threshold, or model outputs can be considered generalisable or clinically applicable.
