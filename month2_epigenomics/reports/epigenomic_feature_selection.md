# Epigenomic Feature Selection Report

## Input dataset

- Patients: 707
- Total input columns: 6,657
- Total methylation features: 6,597
- Non-methylation columns: 60

## Feature-selection method

- Methylation features were converted to numeric values.
- Variance was calculated across patients for every methylation feature.
- Features with missing or zero variance were excluded.
- The top 500 features with the highest variance were retained.

## Output feature table

- Patients: 707
- Final columns: 560
- Selected methylation features: 500

## Top 20 selected methylation features

1. `meth_cg00211215` — variance: 0.109015
2. `meth_cg00103771` — variance: 0.094202
3. `meth_cg00177787` — variance: 0.093156
4. `meth_cg00011616` — variance: 0.090145
5. `meth_cg00084338` — variance: 0.086889
6. `meth_cg00055764` — variance: 0.086230
7. `meth_cg00064255` — variance: 0.085661
8. `meth_cg00063471` — variance: 0.084657
9. `meth_cg00332802` — variance: 0.079924
10. `meth_cg00347563` — variance: 0.078682
11. `meth_cg00100121` — variance: 0.078345
12. `meth_cg00157199` — variance: 0.078049
13. `meth_cg00257271` — variance: 0.077386
14. `meth_cg00197993` — variance: 0.076335
15. `meth_cg00200653` — variance: 0.076104
16. `meth_cg00289546` — variance: 0.075837
17. `meth_cg00204249` — variance: 0.075013
18. `meth_cg00074145` — variance: 0.074742
19. `meth_cg00063174` — variance: 0.074431
20. `meth_cg00065215` — variance: 0.074363

## Output files

- `data/processed/ml_epigenomic_feature_table.csv`
- `outputs/tables/methylation_feature_variance.csv`
- `outputs/tables/selected_methylation_features.csv`
