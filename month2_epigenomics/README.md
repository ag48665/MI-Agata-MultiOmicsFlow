# Month 2 – Epigenomic Integration

## MI-Agata-MultiOmicsFlow

Month 2 extends the clinical-genomic breast cancer analysis developed in Month 1 by integrating DNA methylation data from TCGA-BRCA.

The objective of this stage was to evaluate whether epigenomic information improves predictive performance beyond the clinical and genomic baseline while maintaining leakage-aware model development and validation.

## Workflow

The Month 2 pipeline includes:

- DNA methylation data preparation and quality control
- patient-level matching across clinical, genomic, and epigenomic datasets
- multi-omics feature integration
- methylation feature filtering and selection
- leakage-aware cross-validation
- clinical-genomic versus epigenomic model comparison
- decision-threshold optimization
- probability calibration
- survival modelling
- feature-importance analysis
- final model validation and automated quality control

## Data integration

DNA methylation data were processed and matched to the existing clinical-genomic TCGA-BRCA cohort.

After patient-level matching, the integrated dataset contained:

**707 patients**

with available clinical, genomic, and epigenomic information.

## Predictive performance

Leakage-aware validation indicated substantially improved discrimination after inclusion of methylation-derived features.

| Model | ROC AUC |
|---|---:|
| Clinical + genomic | ~0.548 |
| Clinical + genomic + epigenomic | ~0.866 |

The epigenomic model therefore showed substantially stronger predictive discrimination than the clinical-genomic baseline during internal validation.

## Precision-recall performance

Because the prediction task involved an imbalanced positive class, precision-recall performance was evaluated in addition to ROC AUC.

- **Average precision:** ~0.667
- **Positive-class prevalence:** ~0.124

The average precision was therefore substantially above the prevalence-based reference level.

## Final calibrated model

Following model validation, threshold analysis, and probability calibration, the final development-stage model achieved approximately:

| Metric | Value |
|---|---:|
| ROC AUC | 0.853 |
| PR AUC | 0.660 |
| Brier score | 0.067 |
| Balanced accuracy | 0.779 |
| Recall | 0.773 |
| F1-score | 0.471 |

Isotonic calibration was selected during the calibration stage.

## Threshold optimization

The default classification threshold of 0.50 was evaluated against optimized operating thresholds.

Threshold optimization improved sensitivity and provided a more appropriate operating point for the imbalanced prediction task.

Threshold selection was performed within the validation workflow to reduce the risk of optimistic performance estimates.

## Survival analysis

Survival modelling was performed separately from the binary classification workflow using a Cox proportional-hazards model.

Predicted Cox risk scores were used to stratify patients into lower- and higher-risk groups, and Kaplan-Meier analysis was used to evaluate survival separation between these groups.

## Repository structure

```text
month2_epigenomics/
│
├── scripts/
│   ├── 10_prepare_methylation_data.py
│   ├── 11_clean_methylation_data.py
│   ├── 12_match_methylation_ids.py
│   ├── 13_integrate_epigenomics.py
│   ├── ...
│   └── 29_prepare_submission_materials.py
│
├── figures/
│   └── selected QC, validation, calibration,
│       survival, and model-performance figures
│
├── reports/
│   └── analytical, validation, manuscript,
│       calibration, survival, and QC reports
│
└── tables/
    └── selected model metrics, validation results,
        feature sets, thresholds, and QC summaries
