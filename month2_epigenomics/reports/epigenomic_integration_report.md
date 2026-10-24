# Clinical–Genomic–Epigenomic Integration Report

## Integration summary

- input_methylation_features: 6,597
- input_methylation_samples: 776
- mapped_methylation_samples: 707
- unique_methylation_patients: 707
- unique_clinical_genomic_patients: 1,009
- matched_integrated_patients: 707
- methylation_features_added: 6,597
- duplicate_clinical_patient_rows_removed: 0
- integrated_total_columns: 6,657

## Integrated dataset

- Output file: `clinical_genomic_epigenomic_merged.csv`
- Rows/patients: 707
- Total columns: 6,657

## Processing steps

- Converted the methylation matrix from features × samples to patients × features.
- Converted TCGA sample identifiers to patient identifiers using the Month 2 mapping.
- Averaged methylation values when more than one sample represented the same patient.
- Integrated methylation features with the Month 1 clinical-genomic dataset.
- Retained only patients available in both clinical-genomic and methylation datasets.
