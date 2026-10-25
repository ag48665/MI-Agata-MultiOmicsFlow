# Bulk RNA-seq and Single-cell RNA-seq Integration — Final Summary

## Project
MI-Agata-MultiOmicsFlow — Month 3 Transcriptomics

## Data sources

### Bulk RNA-seq
- Dataset: TCGA-BRCA
- Significant differentially expressed genes: 7,699
- Rows mapped to gene symbols: 6,365
- Unique mapped gene symbols: 6,357

### Single-cell RNA-seq
- Dataset: GSE176078
- Cells analyzed: 9,566
- Biological samples: 26
- Leiden clusters: 25
- Annotated cell-type categories: 9
- Top marker rows: 250
- Unique marker genes: 214

## Integration methodology

1. Loaded significant Bulk RNA-seq genes.
2. Loaded Ensembl-to-gene-symbol mapping.
3. Removed Ensembl version suffixes.
4. Matched Bulk genes to gene symbols.
5. Loaded top 10 Single-cell markers per Leiden cluster.
6. Matched datasets using exact gene symbols.
7. Validated overlap table.
8. Visualized matched genes across Leiden clusters.

## Main results

- Shared unique gene symbols: 80
- Matched cluster-gene rows: 95
- Leiden clusters with matches: 22
- Leiden clusters without matches: 3
- Matched rows with positive Bulk log2FC: 50
- Matched rows with negative Bulk log2FC: 45

The 95 rows represent cluster-gene matches, not 95 independent
genes or biological replicates.

## Quality control

- Exact duplicate rows: 0
- Duplicate cluster-gene pairs: 0
- Missing gene symbols: 0
- Missing Bulk log2FC values: 0
- Missing Bulk adjusted p-values: 0
- Gene-symbol consistency: PASS
- Ambiguous Ensembl mappings detected: 0

## Biological observations

Examples of shared genes include:

- VWF and PECAM1: endothelial-associated markers.
- COL1A1 and COL1A2: extracellular matrix and fibroblast-associated genes.
- MYL9 and TAGLN: contractile cell-associated markers.

These associations are consistent with known marker functions
but do not establish which cell populations caused Bulk RNA-seq
expression differences.

## Limitations

1. The overlap is descriptive and does not constitute
   statistical enrichment analysis.
2. Only the top 10 Single-cell markers per Leiden cluster
   were considered.
3. Not all significant Bulk RNA-seq genes could be mapped
   to gene symbols.
4. Bulk RNA-seq and Single-cell RNA-seq originate from
   different datasets and sample cohorts.
5. Cell-type annotation depends on marker-gene evidence.
6. The analysis does not establish causality or independent
   biological validation.
7. Bulk expression differences can reflect changes in cell
   abundance, cell-intrinsic expression, or both.

## Main outputs

### Tables
- outputs/tables/bulk_single_cell_gene_overlap.csv

### Figures
- outputs/figures/bulk_single_cell_overlap_by_cluster.png

### Reports
- reports/gene_overlap_report.json
- reports/integration_final_summary.md

## Status

Exploratory Bulk RNA-seq and Single-cell RNA-seq gene-overlap
integration completed.

Repository publication review and final reproducibility checks
remain pending.
