# Single-cell RNA-seq — Final Analysis Summary

## Project
MI-Agata-MultiOmicsFlow — Month 3 Transcriptomics

## Dataset
- Accession: GSE176078
- Analysis: Breast cancer single-cell RNA-seq
- Cells after filtering: 9,566
- Biological samples: 26
- Breast cancer subtypes: ER+, HER2+, TNBC

## Analysis workflow
1. Dataset assessment and loading
2. Quality control
3. Cell filtering
4. Normalization
5. Highly variable gene selection and PCA
6. Neighborhood graph, UMAP and Leiden clustering
7. Cluster marker identification
8. Cell-type annotation
9. Cell-type composition by breast cancer subtype
10. Sample-level cell-type composition and QC

## Main results
- 25 Leiden clusters identified.
- 9 annotated cell-type categories.
- Sample distribution:
  - ER+: 11 samples
  - HER2+: 5 samples
  - TNBC: 10 samples
- Sample-level composition table validated.
- All 9,566 analyzed cells accounted for in sample-level composition.

## Main output files
### Figures
- outputs/figures/pca_explained_variance.png
- outputs/figures/umap_leiden.png
- outputs/figures/umap_celltype_annotated.png
- outputs/figures/celltype_composition_by_subtype.png
- outputs/figures/celltype_composition_per_sample.png

### Tables
- outputs/tables/leiden_marker_genes.csv
- outputs/tables/leiden_top10_markers.csv
- outputs/tables/leiden_celltype_annotation_summary.csv
- outputs/tables/celltype_counts_by_subtype.csv
- outputs/tables/celltype_percentages_by_subtype.csv
- outputs/tables/celltype_composition_per_sample.csv

### Quality control
- reports/sample_composition_qc.md

## Interpretation and limitations
Cell-type proportions differ descriptively across breast cancer subtypes.

These observations should not be interpreted as statistically significant
subtype differences without appropriate sample-level statistical testing.

Cell-type annotations are based on marker-gene evidence and require
biological interpretation. Cluster assignments are computational estimates.

The analysis uses 26 biological samples, not 9,566 independent
biological replicates.

## Status
Single-cell RNA-seq exploratory analysis completed.
Final repository security and reproducibility review pending.
