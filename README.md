# Month 3 — Transcriptomics Analysis | TCGA-BRCA

**MI-Agata-MultiOmicsFlow — Transcriptomics Module**

**Dataset:** The Cancer Genome Atlas — Breast Invasive Carcinoma (TCGA-BRCA)  
**Data type:** RNA sequencing (RNA-seq)  
**Analysis:** Differential Gene Expression, Functional Enrichment, Biological Interpretation and Quality Control  
**Status:** Completed — Final QC: PASS

---

## 1. Overview

This module implements a reproducible transcriptomics analysis workflow using RNA-seq data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) cohort.

The workflow includes dataset assessment, RNA-seq data acquisition, gene expression matrix construction, differential gene expression analysis, exploratory visualization, functional enrichment, biological interpretation, and final quality control.

The primary objective is to characterize transcriptional differences in the analyzed TCGA-BRCA comparison and identify biological processes and molecular pathways associated with the observed gene expression changes.

This module forms the transcriptomics component of **MI-Agata-MultiOmicsFlow**, a multi-omics analysis project integrating genomics, epigenomics, transcriptomics, and machine learning approaches for cancer research.

## 2. Research Objectives

The main objectives are to:

1. Acquire, assess, and process TCGA-BRCA RNA-seq data.
2. Construct an expression matrix suitable for downstream analysis.
3. Identify statistically significant differentially expressed genes (DEGs).
4. Explore gene expression patterns using volcano plots, heatmaps, and principal component analysis (PCA).
5. Perform Gene Ontology (GO) Biological Process and KEGG pathway enrichment analysis.
6. Interpret the biological relevance of significant genes and enriched pathways.
7. Validate expected analysis outputs through final quality control.
8. Prepare transcriptomic findings for subsequent multi-omics integration.

## 3. Analysis Workflow

The analysis is organized into sequential Python scripts.

| Step | Script | Description |
|---|---|---|
| 01 | `01_dataset_assessment.py` | Initial dataset assessment |
| 02 | `02_download_tcga_rnaseq.py` | TCGA-BRCA RNA-seq data acquisition |
| 03 | `03_build_expression_matrix.py` | Gene expression matrix construction |
| 04 | `04_differential_expression.py` | Differential gene expression analysis |
| 05 | `05_visualize_differential_expression.py` | Differential expression visualization |
| 06 | `06_heatmap_pca.py` | Heatmap generation and PCA |
| 07 | `07_pathway_enrichment.py` | GO and KEGG functional enrichment |
| 08 | `08_visualize_pathways.py` | Enrichment visualization |
| 09 | `09_biological_interpretation.py` | Biological interpretation |
| 10 | `10_final_qc_summary.py` | Final quality control and summary |

Additional supporting scripts may be used for data preparation and preprocessing.


## 4. Differential Gene Expression Analysis

Differential expression analysis was performed to identify genes exhibiting statistically significant expression changes in the analyzed TCGA-BRCA comparison.

### 4.1 Significance Criteria

Significant differentially expressed genes were selected using:

- Adjusted p-value (FDR) < 0.05
- Absolute log2 fold change > 1

These thresholds combine multiple-testing-adjusted statistical significance with a minimum expression-change criterion.

### 4.2 Differential Expression Results

| Metric | Result |
|---|---:|
| Total genes tested | 28,659 |
| Significant differentially expressed genes | 7,699 |
| Upregulated genes | 2,413 |
| Downregulated genes | 5,286 |

Approximately **26.9% of tested genes** met the differential expression criteria.

Among significant DEGs:

- 31.3% were classified as upregulated.
- 68.7% were classified as downregulated.

The direction of differential expression depends on the comparison and reference groups defined in the analysis.

These results indicate substantial transcriptomic differences within the analyzed comparison and provide the basis for downstream biological interpretation.

### 4.3 Volcano Plot

The volcano plot visualizes gene-level statistical significance against log2 fold change, highlighting genes with large expression changes and strong statistical evidence.

![TCGA-BRCA Volcano Plot](https://raw.githubusercontent.com/ag48665/MI-Agata-MultiOmicsFlow/main/month3_transcriptomics/figures/tcga_brca_volcano_plot.png)



## 5. Exploratory Gene Expression Analysis

Exploratory visualization was used to assess expression patterns, sample relationships, and major sources of variation.

### 5.1 Gene Expression Heatmap

The heatmap visualizes expression patterns among selected genes and samples.

It provides an exploratory overview of expression heterogeneity and potential sample-level patterns.

![TCGA-BRCA Gene Expression Heatmap](https://raw.githubusercontent.com/ag48665/MI-Agata-MultiOmicsFlow/main/month3_transcriptomics/figures/tcga_brca_heatmap.png)



### 5.2 Principal Component Analysis

Principal component analysis (PCA) was used to explore major sources of variation in the expression dataset.

PCA provides a low-dimensional representation of sample relationships and can help identify potential clustering patterns, heterogeneity, or outliers.

![TCGA-BRCA Principal Component Analysis](https://raw.githubusercontent.com/ag48665/MI-Agata-MultiOmicsFlow/main/month3_transcriptomics/figures/tcga_brca_pca.png)


PCA is an exploratory approach and does not independently demonstrate statistically significant group separation.


## 6. Functional Enrichment Analysis

Functional enrichment analysis was performed to identify biological processes and pathways statistically overrepresented among significant differentially expressed genes.

### 6.1 Functional Annotation Databases

Two annotation resources were used:

- **GO Biological Process 2023**
- **KEGG 2021 Human**

### 6.2 Statistical Method

Functional enrichment was performed using a hypergeometric overrepresentation test.

Multiple-testing correction was applied using the Benjamini–Hochberg procedure.

The background gene universe consisted of mapped genes from the differential expression dataset.

This approach evaluates whether genes associated with particular biological annotations are represented more frequently among the selected genes than expected under the statistical null model.

### 6.3 Enrichment Results

| Database | Significant enriched terms |
|---|---:|
| GO Biological Process 2023 | 120 |
| KEGG 2021 Human | 28 |

These results highlight multiple biological processes and molecular pathways associated with the differentially expressed gene set.

## 7. GO Biological Process Enrichment

Prominent enriched biological processes included:

- Extracellular matrix organization
- External encapsulating structure organization
- Extracellular structure organization
- G protein-coupled receptor signaling
- Cell-cell adhesion
- Nervous system development

### 7.1 Biological Interpretation

**Extracellular matrix organization**

The enrichment of extracellular matrix-related processes suggests that genes involved in tissue architecture, extracellular structural organization, and cell–matrix interactions are prominent among the differentially expressed genes.

**Cell adhesion**

Enrichment of cell adhesion-related processes highlights genes involved in cellular interactions and tissue organization.

**Receptor-mediated signaling**

G protein-coupled receptor signaling enrichment indicates representation of genes associated with extracellular signal recognition and cellular communication.

**Developmental processes**

The enrichment of nervous system development-related annotations may reflect shared developmental or signaling genes rather than a direct neurological phenotype.

### 7.2 GO Enrichment Visualization

![GO Biological Process Enrichment](https://raw.githubusercontent.com/ag48665/MI-Agata-MultiOmicsFlow/main/month3_transcriptomics/figures/tcga_brca_go_enrichment.png)

The visualization summarizes selected enriched GO Biological Process terms.

Enrichment results indicate statistical overrepresentation and do not independently establish pathway activation or causal mechanisms.

## 8. KEGG Pathway Enrichment

Significantly enriched KEGG pathways included:

- Neuroactive ligand-receptor interaction
- Cytokine-cytokine receptor interaction
- Complement and coagulation cascades
- Calcium signaling pathway
- ECM-receptor interaction

### 8.1 Biological Interpretation

**Cell signaling**

Neuroactive ligand-receptor interaction and calcium signaling pathways contain genes involved in receptor-mediated communication and intracellular signaling.

**Immune-related processes**

Cytokine-cytokine receptor interaction and complement-related pathways suggest representation of genes associated with immune signaling and inflammatory processes.

**Extracellular matrix interactions**

ECM-receptor interaction enrichment is consistent with the extracellular matrix-related findings identified through GO Biological Process analysis.

Together, these findings identify candidate molecular processes for further investigation in breast cancer research.

### 8.2 KEGG Enrichment Visualization

![KEGG Pathway Enrichment](https://raw.githubusercontent.com/ag48665/MI-Agata-MultiOmicsFlow/main/month3_transcriptomics/figures/tcga_brca_kegg_enrichment.png)

The visualization presents selected enriched KEGG pathways.

Statistical enrichment does not independently demonstrate increased or decreased pathway activity.


## 9. Integrated Biological Interpretation

The differential expression and functional enrichment results highlight several recurring biological themes.

| Biological theme | Supporting findings | Research relevance |
|---|---|---|
| Extracellular matrix organization | GO extracellular matrix processes; KEGG ECM-receptor interaction | Tissue architecture and cell–matrix interactions |
| Receptor-mediated signaling | GPCR signaling; neuroactive ligand-receptor interaction | Cellular communication |
| Calcium signaling | KEGG calcium signaling pathway | Intracellular signal transduction |
| Cell adhesion | GO cell-cell adhesion | Cellular interactions and tissue organization |
| Immune-related signaling | Cytokine-cytokine receptor interaction; complement and coagulation cascades | Immune and inflammatory processes |
| Developmental processes | GO nervous system development | Developmental gene programs represented in the expression dataset |

### Key Findings

The analysis identified widespread differential gene expression and enrichment of biological processes associated with:

1. Extracellular matrix organization and remodeling.
2. Receptor-mediated cellular signaling.
3. Cell adhesion and tissue organization.
4. Immune-related and inflammatory processes.
5. Developmental and signaling-associated gene programs.

These results provide candidate biological themes for further investigation and integration with epigenomic findings.

## 10. Output Files and Repository Structure

The module separates analysis scripts, generated results, visualizations, and reports.

### 10.1 Analysis Scripts

The analysis scripts are located in:

`month3_transcriptomics/scripts/`

The workflow includes dataset assessment, RNA-seq acquisition, expression matrix construction, differential expression, visualization, enrichment, biological interpretation, and quality control.

### 10.2 Figures

The generated visualizations are available in:

`month3_transcriptomics/figures/`

| Figure | Description |
|---|---|
| [Volcano Plot](https://github.com/ag48665/MI-Agata-MultiOmicsFlow/blob/main/month3_transcriptomics/figures/tcga_brca_volcano_plot.png) | Differential expression significance and fold changes |
| [Gene Expression Heatmap](https://github.com/ag48665/MI-Agata-MultiOmicsFlow/blob/main/month3_transcriptomics/figures/tcga_brca_heatmap.png) | Gene expression patterns |
| [PCA Plot](https://github.com/ag48665/MI-Agata-MultiOmicsFlow/blob/main/month3_transcriptomics/figures/tcga_brca_pca.png) | Sample-level variation |
| [GO Enrichment Plot](https://github.com/ag48665/MI-Agata-MultiOmicsFlow/blob/main/month3_transcriptomics/figures/tcga_brca_go_enrichment.png) | GO Biological Process enrichment |
| [KEGG Enrichment Plot](https://github.com/ag48665/MI-Agata-MultiOmicsFlow/blob/main/month3_transcriptomics/figures/tcga_brca_kegg_enrichment.png) | KEGG pathway enrichment |

### 10.3 Main Analysis Tables

Key output tables are located in:

`month3_transcriptomics/outputs/tables/`

Important tables include:

- `tcga_brca_differential_expression.csv`
- `tcga_brca_significant_genes.csv`
- `tcga_brca_top20_genes.csv`
- `tcga_brca_gene_mapping.csv`
- `GO_Biological_Process_2023_enrichment.csv`
- `KEGG_2021_Human_enrichment.csv`
- `tcga_brca_pathway_enrichment.csv`

These tables contain differential expression statistics, gene annotations, functional enrichment results, and summary information.

### 10.4 Reports

Important reports include:

- `month3_transcriptomics/reports/tcga_brca_biological_interpretation.md`
- `month3_transcriptomics/reports/tcga_brca_final_summary.md`
- `month3_transcriptomics/reports/tcga_brca_final_qc.csv`

These reports document the biological interpretation, overall analysis summary, and final quality control results.

## 11. Quality Control

Final quality control was performed to assess the presence and validity of expected analysis outputs according to the implemented validation checks.

### 11.1 QC Summary

| Quality control metric | Result |
|---|---:|
| Expected files checked | 16 |
| Passed checks | 16 |
| Failed checks | 0 |
| Final QC status | PASS |

All 16 expected output files passed the implemented quality control checks.

### 11.2 Interpretation of QC Results

The final QC result supports the completeness of the expected analysis outputs under the implemented validation criteria.

However, successful file-level QC does not independently establish:

- Biological validity of all findings.
- Correctness of every statistical assumption.
- Absence of confounding factors.
- External reproducibility across computational environments.
- Generalizability to independent patient cohorts.

These considerations remain important when interpreting the results and planning downstream analyses.

## 12. Software and Dependencies

The workflow uses Python and scientific computing libraries, including:

| Library | Main purpose |
|---|---|
| `pandas` | Data processing and tabular analysis |
| `numpy` | Numerical computation |
| `scipy` | Statistical analysis |
| `matplotlib` | Scientific visualization |
| `statsmodels` | Statistical modeling and multiple-testing procedures |
| `gseapy` | Functional enrichment analysis |

A dedicated Python virtual environment is used for the transcriptomics module.

### Reproducibility Considerations

Reproducibility depends on maintaining consistent:

- Input datasets and data versions.
- Preprocessing procedures.
- Gene annotation and mapping methods.
- Differential expression methodology.
- Statistical thresholds.
- Functional annotation database versions.
- Software dependencies and computational environment.

The repository provides analysis scripts and generated outputs to support workflow transparency and reproducibility.

Raw or restricted datasets may be excluded from the public repository for data governance, licensing, privacy, or storage reasons.

## 13. Scientific Limitations

Several considerations are important when interpreting the results.

### 13.1 Cohort Composition

Observed expression differences may be influenced by sample composition, tumor heterogeneity, clinical characteristics, and other potential confounding factors.

### 13.2 Statistical Interpretation

Differential expression identifies statistical associations between expression measurements and the analyzed comparison.

It does not independently establish causal relationships.

### 13.3 Functional Enrichment

Enrichment analysis is sensitive to gene selection criteria, annotation coverage, background gene definitions, and pathway database structure.

Overlapping gene sets may generate redundant enrichment findings.

### 13.4 Pathway Activity

Overrepresentation of genes within a pathway does not directly demonstrate pathway activation or suppression.

Additional analyses are required to investigate directional pathway changes.

### 13.5 External Validation

Independent datasets and additional experimental or computational validation may be required to assess the robustness and generalizability of the findings.

## 14. Future Multi-Omics Integration

The next stage of MI-Agata-MultiOmicsFlow will focus on integrating transcriptomic findings with genomic and epigenomic results.

Planned analyses include:

1. Identification of genes supported by multiple omics layers.
2. Comparison of differential gene expression and DNA methylation findings.
3. Exploration of potential relationships between epigenomic changes and gene expression.
4. Integrated functional enrichment and pathway interpretation.
5. Visualization of cross-omics molecular relationships.
6. Assessment of potential multi-omics features for downstream predictive modeling.
7. Extension of reproducibility and project-level quality control documentation.

The objective is to develop a more comprehensive molecular interpretation than can be obtained from a single omics layer alone.

## 15. Conclusion

The Month 3 transcriptomics module successfully completed the implemented TCGA-BRCA RNA-seq analysis workflow and final output quality control.

The analysis identified **7,699 significant differentially expressed genes among 28,659 tested genes**, together with **120 significant GO Biological Process terms and 28 significant KEGG pathways**.

The findings highlight extracellular matrix organization, receptor-mediated signaling, cell adhesion, and immune-related biological processes as prominent themes within the analyzed differential expression results.

These observations provide a foundation for hypothesis generation and subsequent integration with epigenomic and genomic findings.

The module contributes a structured transcriptomics analysis component to the broader MI-Agata-MultiOmicsFlow research pipeline.

---

**Project:** [MI-Agata-MultiOmicsFlow](https://github.com/ag48665/MI-Agata-MultiOmicsFlow)  
**Module:** Month 3 — Transcriptomics  
**Dataset:** TCGA-BRCA  
**Final analysis status:** Completed  
**Final implemented QC status:** PASS
