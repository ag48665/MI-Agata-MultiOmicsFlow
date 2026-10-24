# TCGA-BRCA Biological Interpretation

## 1. Study overview

RNA-seq differential expression analysis comparing primary breast tumors with normal breast tissue.

Differential expression was evaluated using DESeq2. Pathway over-representation analysis used GO Biological Process and KEGG gene sets.

Significance thresholds: adjusted p-value < 0.05 and absolute log2 fold change > 1 for gene selection.

## 2. Differential expression results

- Genes tested: 28659
- Significant genes: 7699
- Upregulated in tumor: 2413
- Downregulated in tumor: 5286

### Top 20 genes by adjusted p-value

```
           gene_id log2FoldChange       padj
ENSG00000099953.10      6.110e+00  0.000e+00
ENSG00000123500.10      7.571e+00  0.000e+00
 ENSG00000230838.1      6.547e+00 6.121e-274
ENSG00000119771.15     -2.981e+00 7.098e-236
ENSG00000060718.22      6.098e+00 2.660e-233
 ENSG00000165197.5     -5.729e+00 1.101e-229
ENSG00000137745.12      7.304e+00 7.306e-212
ENSG00000203805.11      5.564e+00 8.142e-212
ENSG00000117650.13      4.549e+00 1.376e-202
 ENSG00000154736.6     -3.131e+00 1.707e-200
ENSG00000137225.13     -3.052e+00 2.593e-200
ENSG00000122641.11      3.713e+00 1.575e-198
ENSG00000090889.12      4.053e+00 2.980e-197
ENSG00000249669.10     -3.895e+00 1.069e-196
ENSG00000167434.10     -7.111e+00 9.863e-196
 ENSG00000286289.1     -3.192e+00 2.188e-190
 ENSG00000241684.6     -3.589e+00 1.532e-189
ENSG00000119147.10     -5.057e+00 1.093e-182
ENSG00000155886.11      5.840e+00 1.534e-175
ENSG00000136158.12     -2.645e+00 2.098e-175
```

## 3. GO Biological Process enrichment

Significant GO terms: 120

### Top 10 GO terms

```
                                                                                  Term  Overlap Adjusted P-value
                                        Extracellular Matrix Organization (GO:0030198)       94        3.531e-12
                            External Encapsulating Structure Organization (GO:0045229)       64        2.331e-10
                                     Extracellular Structure Organization (GO:0043062)       63        6.748e-10
Adenylate Cyclase-Modulating G Protein-Coupled Receptor Signaling Pathway (GO:0007188)       76        6.531e-08
                                           Antimicrobial Humoral Response (GO:0019730)       45        4.889e-07
                Cell-Cell Adhesion Via Plasma-Membrane Adhesion Molecules (GO:0098742)       79        1.555e-06
                                                    Epidermis Development (GO:0008544)       46        1.555e-06
                                                      Metal Ion Transport (GO:0030001)       77        1.555e-06
Adenylate Cyclase-Activating G Protein-Coupled Receptor Signaling Pathway (GO:0007189)       52        6.312e-06
                                               Nervous System Development (GO:0007399)      161        8.380e-06
```

## 4. KEGG pathway enrichment

Significant KEGG pathways: 28

### Top 10 KEGG pathways

```
                                                         Term  Overlap Adjusted P-value
                      Neuroactive ligand-receptor interaction      136        5.310e-14
                                 Systemic lupus erythematosus       57        1.392e-08
                       Cytokine-cytokine receptor interaction      109        2.550e-07
                                                   Alcoholism       71        2.958e-06
                          Complement and coagulation cascades       41        4.216e-06
                                           Nicotine addiction       24        1.113e-05
                                           Retinol metabolism       33        1.236e-05
Viral protein interaction with cytokine and cytokine receptor       47        3.196e-05
                                    Calcium signaling pathway       90        1.468e-04
                      Neutrophil extracellular trap formation       66        4.237e-04
```

## 5. Biological interpretation

The enriched gene sets highlight biological processes that are over-represented among differentially expressed genes in the tumor-versus-normal comparison.

Potential themes include extracellular matrix remodeling, cell adhesion, immune signaling, and receptor-mediated communication. These themes should be interpreted in the context of the actual enriched terms.

Enrichment does not establish pathway activation or inhibition, and does not demonstrate causality.

KEGG categories such as Alcoholism or Nicotine addiction refer to annotated gene sets and do not imply patient substance use.

## 6. Limitations

- Tumor and normal samples may be paired by patient; the statistical design should account for pairing.
- Results depend on gene identifier mapping and the reference gene-set libraries.
- Over-representation analysis does not measure directional pathway activity.
- Findings require independent biological validation.

## 7. Conclusion

The analysis identifies statistically significant transcriptional differences between TCGA-BRCA tumor and normal samples and highlights candidate biological processes for further investigation.