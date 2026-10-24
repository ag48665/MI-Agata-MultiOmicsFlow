# Epigenomic Survival Model Report

## Cohort

- Total patients: 706
- Training patients: 564
- Testing patients: 142
- Events/deaths: 88
- Censored/living: 618

## Survival outcome

- Duration variable: `os_months`
- Event variable: `outcome`
- Event = 1: deceased
- Event = 0: censored/living

## Predictors

- Final encoded predictors: 170
- Selected methylation features: 50

## Model

- Algorithm: Cox proportional hazards regression
- Penalizer: 1.0

## Performance

- Training C-index: 0.917
- Testing C-index: 0.904

## Interpretation

A concordance index close to 0.5 indicates approximately random risk ranking. Higher values indicate better discrimination between patients with shorter and longer survival.

Hazard ratios greater than 1 indicate increased estimated hazard, whereas hazard ratios below 1 indicate decreased estimated hazard.

## Leakage controls

- Patient and sample identifiers were excluded from model predictors.
- Other survival-related variables were excluded from predictors.
- High-cardinality categorical variables were removed.
- Methylation feature selection was performed using training data only.
- Missing-value imputation was based on training data only.
