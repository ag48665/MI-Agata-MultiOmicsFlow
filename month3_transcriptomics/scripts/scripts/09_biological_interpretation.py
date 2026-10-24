
from pathlib import Path
import pandas as pd
import numpy as np

# ==========================================
# 09 — BIOLOGICAL INTERPRETATION
# TCGA-BRCA Transcriptomics
# ==========================================

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "outputs" / "tables"
REPORTS = ROOT / "reports"

REPORTS.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

DE_FILE = TABLES / "tcga_brca_differential_expression.csv"
GO_FILE = TABLES / "GO_Biological_Process_2023_enrichment.csv"
KEGG_FILE = TABLES / "KEGG_2021_Human_enrichment.csv"

FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1


def load_results():
    print("=== LOADING RESULTS ===")

    for path in [DE_FILE, GO_FILE, KEGG_FILE]:
        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

    de = pd.read_csv(DE_FILE)
    go = pd.read_csv(GO_FILE)
    kegg = pd.read_csv(KEGG_FILE)

    required_de = {"gene_id", "padj", "log2FoldChange"}
    required_pathway = {"Term", "Adjusted P-value", "Overlap"}

    if not required_de.issubset(de.columns):
        raise ValueError("Missing DESeq2 columns")

    for name, data in [("GO", go), ("KEGG", kegg)]:
        if not required_pathway.issubset(data.columns):
            raise ValueError(f"Missing {name} columns")

    return de, go, kegg


def analyze_genes(de):
    print("\n=== DIFFERENTIAL EXPRESSION ===")

    de = de.copy()
    de["padj"] = pd.to_numeric(de["padj"], errors="coerce")
    de["log2FoldChange"] = pd.to_numeric(
        de["log2FoldChange"], errors="coerce"
    )

    de = de.dropna(subset=["padj", "log2FoldChange"])

    significant = de[
        (de["padj"] < FDR_THRESHOLD)
        & (de["log2FoldChange"].abs() > LOG2FC_THRESHOLD)
    ].copy()

    up = significant[significant["log2FoldChange"] > 0]
    down = significant[significant["log2FoldChange"] < 0]

    top_genes = significant.sort_values(
        ["padj", "gene_id"], ascending=[True, True]
    ).head(20)

    output = TABLES / "tcga_brca_top20_genes.csv"
    top_genes.to_csv(output, index=False)

    summary = {
        "Genes tested": len(de),
        "Significant genes": len(significant),
        "Upregulated in tumor": len(up),
        "Downregulated in tumor": len(down),
    }

    for key, value in summary.items():
        print(f"{key}: {value}")

    print("Saved:", output)

    return summary, top_genes


def analyze_pathways(df, library):
    print(f"\n=== {library} ===")

    df = df.copy()

    df["Adjusted P-value"] = pd.to_numeric(
        df["Adjusted P-value"], errors="coerce"
    )
    df["Overlap"] = pd.to_numeric(
        df["Overlap"], errors="coerce"
    )

    df = df.dropna(
        subset=["Term", "Adjusted P-value", "Overlap"]
    )

    significant = df[
        (df["Adjusted P-value"] >= 0)
        & (df["Adjusted P-value"] < FDR_THRESHOLD)
    ].copy()

    significant = significant.sort_values(
        ["Adjusted P-value", "Term"]
    )

    top = significant.head(10)

    output = TABLES / f"tcga_brca_{library.lower()}_top10.csv"
    top.to_csv(output, index=False)

    print("Significant pathways:", len(significant))
    print("Saved:", output)

    return len(significant), top


def format_table(df, columns):
    if df.empty:
        return "No significant results."

    selected = df[columns].copy()

    for col in selected.columns:
        if pd.api.types.is_float_dtype(selected[col]):
            selected[col] = selected[col].map(
                lambda x: f"{x:.3e}" if pd.notna(x) else "NA"
            )

    return selected.to_string(index=False)


def generate_report(gene_summary, top_genes,
                    go_count, go_top,
                    kegg_count, kegg_top):

    print("\n=== GENERATING REPORT ===")

    report = []

    report.append("# TCGA-BRCA Biological Interpretation")
    report.append("")
    report.append("## 1. Study overview")
    report.append("")
    report.append(
        "RNA-seq differential expression analysis comparing "
        "primary breast tumors with normal breast tissue."
    )
    report.append("")
    report.append(
        "Differential expression was evaluated using DESeq2. "
        "Pathway over-representation analysis used GO Biological "
        "Process and KEGG gene sets."
    )
    report.append("")
    report.append(
        "Significance thresholds: adjusted p-value < 0.05 "
        "and absolute log2 fold change > 1 for gene selection."
    )

    report.append("")
    report.append("## 2. Differential expression results")
    report.append("")

    for key, value in gene_summary.items():
        report.append(f"- {key}: {value}")

    report.append("")
    report.append("### Top 20 genes by adjusted p-value")
    report.append("")
    report.append("```")
    report.append(
        format_table(
            top_genes,
            ["gene_id", "log2FoldChange", "padj"]
        )
    )
    report.append("```")

    report.append("")
    report.append("## 3. GO Biological Process enrichment")
    report.append("")
    report.append(f"Significant GO terms: {go_count}")
    report.append("")
    report.append("### Top 10 GO terms")
    report.append("")
    report.append("```")
    report.append(
        format_table(
            go_top,
            ["Term", "Overlap", "Adjusted P-value"]
        )
    )
    report.append("```")

    report.append("")
    report.append("## 4. KEGG pathway enrichment")
    report.append("")
    report.append(f"Significant KEGG pathways: {kegg_count}")
    report.append("")
    report.append("### Top 10 KEGG pathways")
    report.append("")
    report.append("```")
    report.append(
        format_table(
            kegg_top,
            ["Term", "Overlap", "Adjusted P-value"]
        )
    )
    report.append("```")

    report.append("")
    report.append("## 5. Biological interpretation")
    report.append("")
    report.append(
        "The enriched gene sets highlight biological processes "
        "that are over-represented among differentially expressed "
        "genes in the tumor-versus-normal comparison."
    )
    report.append("")
    report.append(
        "Potential themes include extracellular matrix remodeling, "
        "cell adhesion, immune signaling, and receptor-mediated "
        "communication. These themes should be interpreted "
        "in the context of the actual enriched terms."
    )
    report.append("")
    report.append(
        "Enrichment does not establish pathway activation or "
        "inhibition, and does not demonstrate causality."
    )
    report.append("")
    report.append(
        "KEGG categories such as Alcoholism or Nicotine addiction "
        "refer to annotated gene sets and do not imply patient "
        "substance use."
    )

    report.append("")
    report.append("## 6. Limitations")
    report.append("")
    report.append(
        "- Tumor and normal samples may be paired by patient; "
        "the statistical design should account for pairing."
    )
    report.append(
        "- Results depend on gene identifier mapping and "
        "the reference gene-set libraries."
    )
    report.append(
        "- Over-representation analysis does not measure "
        "directional pathway activity."
    )
    report.append(
        "- Findings require independent biological validation."
    )

    report.append("")
    report.append("## 7. Conclusion")
    report.append("")
    report.append(
        "The analysis identifies statistically significant "
        "transcriptional differences between TCGA-BRCA tumor "
        "and normal samples and highlights candidate biological "
        "processes for further investigation."
    )

    output = REPORTS / "tcga_brca_biological_interpretation.md"
    output.write_text("\n".join(report), encoding="utf-8")

    print("Saved:", output)


def main():
    print("=" * 55)
    print("TCGA-BRCA BIOLOGICAL INTERPRETATION")
    print("=" * 55)

    de, go, kegg = load_results()

    gene_summary, top_genes = analyze_genes(de)

    go_count, go_top = analyze_pathways(
        go, "GO_Biological_Process_2023"
    )

    kegg_count, kegg_top = analyze_pathways(
        kegg, "KEGG_2021_Human"
    )

    generate_report(
        gene_summary,
        top_genes,
        go_count,
        go_top,
        kegg_count,
        kegg_top
    )

    print("\nSTATUS: PASS")
    print("BIOLOGICAL INTERPRETATION COMPLETED")


if __name__ == "__main__":
    main()
