# Single-cell RNA-seq — Sample Composition QC

## Dataset
- GEO accession: GSE176078
- Cells after filtering: 9,566
- Biological samples: 26
- Cell types: 9
- Breast cancer subtypes:
  - ER+: 11 samples
  - HER2+: 5 samples
  - TNBC: 10 samples

## Validation results

| QC check | Result |
|---|---|
| 26 unique samples | PASS |
| 9 cell types per sample | PASS |
| Cell counts match sample totals | PASS |
| Percentages sum to 100% | PASS |
| No negative cell counts | PASS |
| No duplicate sample-celltype pairs | PASS |
| Total cells = 9,566 | PASS |

Overall result: 7 PASS / 0 FAIL.

## Interpretation
The sample-level cell composition table passed all seven
structural and arithmetic validation checks.

The analysis uses biological samples as the appropriate
unit for subsequent comparisons between breast cancer
subtypes.

Differences between ER+, HER2+, and TNBC should not be
interpreted as statistically significant without
appropriate sample-level statistical testing.

## Output files
- outputs/tables/celltype_composition_per_sample.csv
- outputs/figures/celltype_composition_per_sample.png
- reports/sample_composition_report.json
