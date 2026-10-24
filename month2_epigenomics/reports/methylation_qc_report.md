# TCGA-BRCA Methylation QC Report

## Dataset dimensions

- Input rows/features: 6,597
- Cleaned rows/features: 6,597
- Total columns: 780
- Sample columns: 776
- Metadata columns: 4

## Metadata columns

```text
ENTITY_STABLE_ID
NAME
DESCRIPTION
TRANSCRIPT_ID
```

## Quality-control results

- Total methylation values: 5,119,272
- Missing values: 14,070
- Missing percentage: 0.2748%
- Duplicated rows: 0
- Values below 0: 0
- Values above 1: 0
- Minimum value: 0.0
- Maximum value: 0.995715

## Cleaning actions

- Converted sample columns to numeric values.
- Converted invalid entries to missing values.
- Removed fully duplicated rows.
- Restricted methylation values to the 0–1 range.
- Retained missing values for later imputation.
