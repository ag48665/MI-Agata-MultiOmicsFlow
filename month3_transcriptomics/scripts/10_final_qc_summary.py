
from pathlib import Path
import pandas as pd
import numpy as np

# ==========================================
# 10 - FINAL QUALITY CONTROL & SUMMARY
# TCGA-BRCA | Month 3 Transcriptomics
# ==========================================

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "figures"
REPORTS = ROOT / "reports"

REPORTS.mkdir(parents=True, exist_ok=True)

REQUIRED_TABLES = [
    "tcga_brca_differential_expression.csv",
    "tcga_brca_significant_genes.csv",
    "tcga_brca_pca_coordinates.csv",
    "tcga_brca_gene_mapping.csv",
    "GO_Biological_Process_2023_enrichment.csv",
    "KEGG_2021_Human_enrichment.csv",
    "tcga_brca_pathway_enrichment.csv",
    "tcga_brca_top20_genes.csv",
    "tcga_brca_go_biological_process_2023_top10.csv",
    "tcga_brca_kegg_2021_human_top10.csv",
]

REQUIRED_FIGURES = [
    "tcga_brca_volcano_plot.png",
    "tcga_brca_heatmap.png",
    "tcga_brca_pca.png",
    "tcga_brca_go_enrichment.png",
    "tcga_brca_kegg_enrichment.png",
]

REQUIRED_REPORTS = [
    "tcga_brca_biological_interpretation.md",
]


def check_files(folder, filenames):
    results = []

    for filename in filenames:
        path = folder / filename
        exists = path.is_file()
        size = path.stat().st_size if exists else 0

        status = "PASS" if exists and size > 0 else "FAIL"

        results.append({
            "File": filename,
            "Status": status,
            "Size_bytes": size
        })

        print(f"[{status}] {filename}")

    return pd.DataFrame(results)


def main():
    print("=" * 55)
    print("TCGA-BRCA FINAL QC & SUMMARY")
    print("=" * 55)

    print("\n=== CHECKING TABLES ===")
    tables_qc = check_files(TABLES, REQUIRED_TABLES)

    print("\n=== CHECKING FIGURES ===")
    figures_qc = check_files(FIGURES, REQUIRED_FIGURES)

    print("\n=== CHECKING REPORTS ===")
    reports_qc = check_files(REPORTS, REQUIRED_REPORTS)

    all_qc = pd.concat(
        [tables_qc, figures_qc, reports_qc],
        ignore_index=True
    )

    qc_output = REPORTS / "tcga_brca_final_qc.csv"
    all_qc.to_csv(qc_output, index=False)

    summary = [
        "# TCGA-BRCA | Month 3 Final QC Summary",
        "",
        "## File validation",
        f"- Files checked: {len(all_qc)}",
        f"- Passed: {(all_qc['Status'] == 'PASS').sum()}",
        f"- Failed: {(all_qc['Status'] == 'FAIL').sum()}",
    ]

    de_file = TABLES / "tcga_brca_differential_expression.csv"

    if de_file.is_file():
        de = pd.read_csv(de_file)

        required = {"padj", "log2FoldChange"}

        if required.issubset(de.columns):
            de["padj"] = pd.to_numeric(
                de["padj"], errors="coerce"
            )
            de["log2FoldChange"] = pd.to_numeric(
                de["log2FoldChange"], errors="coerce"
            )

            valid = de.dropna(
                subset=["padj", "log2FoldChange"]
            )

            sig = valid[
                (valid["padj"] < 0.05)
                & (valid["log2FoldChange"].abs() > 1)
            ]

            up = (sig["log2FoldChange"] > 0).sum()
            down = (sig["log2FoldChange"] < 0).sum()

            summary += [
                "",
                "## Differential expression",
                f"- Genes in results: {len(de)}",
                f"- Genes with valid statistics: {len(valid)}",
                f"- Significant genes: {len(sig)}",
                f"- Upregulated: {up}",
                f"- Downregulated: {down}",
            ]

            print("\n=== DIFFERENTIAL EXPRESSION ===")
            print("Genes tested:", len(de))
            print("Significant genes:", len(sig))
            print("Upregulated:", up)
            print("Downregulated:", down)
        else:
            summary += [
                "",
                "## Differential expression",
                "- Required statistical columns missing."
            ]

    for library in [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
    ]:
        path = TABLES / f"{library}_enrichment.csv"

        if path.is_file():
            data = pd.read_csv(path)

            if "Adjusted P-value" in data.columns:
                fdr = pd.to_numeric(
                    data["Adjusted P-value"],
                    errors="coerce"
                )
                count = (fdr < 0.05).sum()

                summary += [
                    "",
                    f"## {library}",
                    f"- Significant pathways: {count}",
                ]

                print(f"{library}: {count} significant pathways")

    failures = (all_qc["Status"] == "FAIL").sum()

    summary += [
        "",
        "## Final file QC status",
        "PASS" if failures == 0 else "INCOMPLETE",
        "",
        "Note: File QC checks presence and non-empty files.",
        "It does not independently validate biological conclusions.",
    ]

    report_path = REPORTS / "tcga_brca_final_summary.md"
    report_path.write_text(
        "\n".join(summary),
        encoding="utf-8"
    )

    print("\n=== FINAL SUMMARY ===")
    print("Files checked:", len(all_qc))
    print("Passed:", (all_qc["Status"] == "PASS").sum())
    print("Failed:", failures)
    print("QC table:", qc_output)
    print("Summary:", report_path)

    if failures == 0:
        print("\nSTATUS: PASS")
    else:
        print("\nSTATUS: INCOMPLETE")


if __name__ == "__main__":
    main()
