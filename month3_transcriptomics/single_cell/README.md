# Single-cell RNA-seq — Breast Cancer

## Project

MI-Agata-MultiOmicsFlow — Month 3: Transcriptomics

## Overview

This module contains an exploratory single-cell RNA-seq
analysis of breast cancer using the GEO GSE176078 dataset.

The workflow includes:
- Dataset assessment
- Quality control and filtering
- Normalization
- PCA and UMAP
- Leiden clustering
- Marker gene analysis
- Cell-type annotation
- Cell-type composition analysis

## Results

The main analysis summary is available in:

reports/single_cell_final_summary.md

## Environment

Python 3.13.14
Scanpy 1.12.4

Dependencies are documented in:
requirements-single-cell.txt

## Integration

Exploratory integration with Bulk RNA-seq results is
documented in ../integration/.

The analysis identified 80 shared unique gene symbols
across 95 cluster-gene overlap rows and 22 Leiden clusters.

The overlap is descriptive and does not constitute
independent biological validation.

## Limitations

The analysis is exploratory. Cell-type annotations
and gene overlap results require cautious interpretation.