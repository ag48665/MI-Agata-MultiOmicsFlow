# Model Validation and Data Leakage Audit

## Objective

This validation stage was created to investigate the unusually high performance observed in the initial Month 2 epigenomic classifier.

The validation uses paired stratified cross-validation and performs methylation feature selection independently inside each training fold.

## Cohort

- Matched patients: 707
- Outcome-positive patients: 88
- Outcome-negative patients: 619

## Cross-validation design

- Stratified folds: 5
- Random Forest trees per fold: 300
- Identical patient folds were used for clinical-genomic and epigenomic models.
- Missing-value imputation was fitted using training data only.
- Categorical encoding was fitted using training data only.
- High-cardinality categorical variables were excluded.
- Methylation variance selection was performed using training data only.

## Feature sets

- Clinical-genomic candidate features: 17
- Non-methylation Month 2 candidate features: 47
- Available methylation features: 6,597
- Methylation features selected per fold: 500

## Leakage audit

- Features excluded as identifiers, survival variables or target proxies: 16

| dataset                     | feature                 | reason                        |
|:----------------------------|:------------------------|:------------------------------|
| clinical_genomic            | patient_id              | exact_exclusion               |
| clinical_genomic            | outcome                 | exact_exclusion               |
| clinical_genomic_epigenomic | patient_id              | exact_exclusion               |
| clinical_genomic_epigenomic | other_patient_id        | identifier:patient_id         |
| clinical_genomic_epigenomic | os_status               | exact_exclusion               |
| clinical_genomic_epigenomic | os_months               | survival_time                 |
| clinical_genomic_epigenomic | dss_status              | dss_variable                  |
| clinical_genomic_epigenomic | dss_months              | survival_time                 |
| clinical_genomic_epigenomic | dfs_status              | dfs_variable                  |
| clinical_genomic_epigenomic | dfs_months              | survival_time                 |
| clinical_genomic_epigenomic | pfs_status              | pfs_variable                  |
| clinical_genomic_epigenomic | pfs_months              | survival_time                 |
| clinical_genomic_epigenomic | sample_id               | exact_exclusion               |
| clinical_genomic_epigenomic | tissue_source_site_code | identifier:tissue_source_site |
| clinical_genomic_epigenomic | tissue_source_site      | identifier:tissue_source_site |
| clinical_genomic_epigenomic | outcome                 | exact_exclusion               |

## Cross-validation summary

| model                       |   folds |   accuracy_mean |   accuracy_std |   balanced_accuracy_mean |   balanced_accuracy_std |   precision_mean |   precision_std |   recall_mean |   recall_std |   f1_score_mean |   f1_score_std |   roc_auc_mean |   roc_auc_std |   roc_auc_change_vs_baseline |
|:----------------------------|--------:|----------------:|---------------:|-------------------------:|------------------------:|-----------------:|----------------:|--------------:|-------------:|----------------:|---------------:|---------------:|--------------:|-----------------------------:|
| clinical_genomic            |       5 |          0.7963 |         0.0388 |                   0.4794 |                  0.0384 |           0.0954 |          0.0926 |        0.0575 |       0.0416 |          0.0708 |         0.0573 |         0.5482 |        0.0264 |                       0      |
| clinical_genomic_epigenomic |       5 |          0.8953 |         0.0095 |                   0.5835 |                  0.0396 |           0.9    |          0.2236 |        0.1686 |       0.0763 |          0.2813 |         0.1187 |         0.8657 |        0.0505 |                       0.3175 |

## Fold-level results

| model                       |   fold |   accuracy |   balanced_accuracy |   f1_score |   roc_auc |
|:----------------------------|-------:|-----------:|--------------------:|-----------:|----------:|
| clinical_genomic            |      1 |     0.7817 |              0.4713 |     0.0606 |    0.5461 |
| clinical_genomic_epigenomic |      1 |     0.9014 |              0.6111 |     0.3636 |    0.9444 |
| clinical_genomic            |      2 |     0.7958 |              0.4794 |     0.0645 |    0.5513 |
| clinical_genomic_epigenomic |      2 |     0.9014 |              0.6111 |     0.3636 |    0.8058 |
| clinical_genomic            |      3 |     0.8511 |              0.5346 |     0.16   |    0.5906 |
| clinical_genomic_epigenomic |      3 |     0.8936 |              0.5588 |     0.2105 |    0.8712 |
| clinical_genomic            |      4 |     0.8085 |              0.4851 |     0.069  |    0.5218 |
| clinical_genomic_epigenomic |      4 |     0.8794 |              0.5254 |     0.1053 |    0.8468 |
| clinical_genomic            |      5 |     0.7447 |              0.4268 |     0      |    0.5312 |
| clinical_genomic_epigenomic |      5 |     0.9007 |              0.6111 |     0.3636 |    0.8604 |

## Methylation feature-selection stability

The following table shows features most frequently selected among the five independent training folds.

| feature         |   folds_selected |   selection_frequency |
|:----------------|-----------------:|----------------------:|
| meth_cg00163018 |                5 |                     1 |
| meth_cg00254681 |                5 |                     1 |
| meth_cg00321288 |                5 |                     1 |
| meth_cg00268185 |                5 |                     1 |
| meth_cg00217080 |                5 |                     1 |
| meth_cg00244745 |                5 |                     1 |
| meth_cg00157796 |                5 |                     1 |
| meth_cg00343633 |                5 |                     1 |
| meth_cg00377087 |                5 |                     1 |
| meth_cg00002719 |                5 |                     1 |
| meth_cg00204249 |                5 |                     1 |
| meth_cg00241712 |                5 |                     1 |
| meth_cg00199007 |                5 |                     1 |
| meth_cg00347904 |                5 |                     1 |
| meth_cg00197993 |                5 |                     1 |
| meth_cg00065215 |                5 |                     1 |
| meth_cg00063174 |                5 |                     1 |
| meth_cg00268840 |                5 |                     1 |
| meth_cg00074145 |                5 |                     1 |
| meth_cg00376544 |                5 |                     1 |

## Interpretation

If the cross-validated epigenomic ROC AUC remains close to 1.0 across independent folds after the leakage controls above, the strong predictive signal is more credible.

If performance decreases substantially relative to the original hold-out result, the original ROC AUC of 1.0 should not be reported as a stable estimate of generalisation performance.

Balanced accuracy, recall and F1-score should also be considered because the outcome distribution is imbalanced.
