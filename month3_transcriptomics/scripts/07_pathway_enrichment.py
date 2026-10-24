
from pathlib import Path
import pandas as pd
import gseapy as gp
import mygene

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1. PATHS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

TABLES = ROOT / "outputs" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

INPUT = TABLES / "tcga_brca_differential_expression.csv"

if not INPUT.exists():
    raise FileNotFoundError(f"Missing input file: {INPUT}")


# ============================================================
# 2. LOAD DESEQ2 RESULTS
# ============================================================

print("\n=== LOADING DESEQ2 RESULTS ===")

df = pd.read_csv(INPUT)

required = {"gene_id", "padj", "log2FoldChange"}

if not required.issubset(df.columns):
    raise ValueError(
        f"Missing columns: {required - set(df.columns)}"
    )

df = df.dropna(
    subset=["gene_id", "padj", "log2FoldChange"]
).copy()

df["ensembl_id"] = (
    df["gene_id"].astype(str).str.split(".").str[0]
)

df = df.drop_duplicates("ensembl_id")

significant = df[
    (df["padj"] < 0.05) &
    (df["log2FoldChange"].abs() > 1)
].copy()

print("Genes tested:", len(df))
print("Significant genes:", len(significant))


# ============================================================
# 3. MAP ENSEMBL IDS TO GENE SYMBOLS
# ============================================================

print("\n=== GENE ID MAPPING ===")

mg = mygene.MyGeneInfo()

gene_ids = df["ensembl_id"].dropna().unique().tolist()

print("Mapping genes using MyGene.info...")

mapping_results = mg.querymany(
    gene_ids,
    scopes="ensembl.gene",
    fields="symbol",
    species="human",
    as_dataframe=False,
    verbose=False
)

mapping = pd.DataFrame(mapping_results)

if not {"query", "symbol"}.issubset(mapping.columns):
    raise ValueError(
        "MyGene.info returned no usable gene mapping."
    )

mapping = mapping.rename(
    columns={"query": "ensembl_id"}
)

mapping = mapping.dropna(
    subset=["ensembl_id", "symbol"]
)

mapping = mapping[
    ["ensembl_id", "symbol"]
].copy()

mapping["symbol"] = mapping["symbol"].astype(str)

mapping = mapping[
    mapping["symbol"].str.strip() != ""
]

mapping = mapping.drop_duplicates("ensembl_id")

print("Mapped Ensembl IDs:", len(mapping))

mapping.to_csv(
    TABLES / "tcga_brca_gene_mapping.csv",
    index=False
)


# ============================================================
# 4. PREPARE BACKGROUND AND SIGNIFICANT GENES
# ============================================================

print("\n=== PREPARING GENE SETS ===")

background = df.merge(
    mapping,
    on="ensembl_id",
    how="inner"
)

foreground = significant.merge(
    mapping,
    on="ensembl_id",
    how="inner"
)

background_genes = set(background["symbol"])

foreground_genes = (
    set(foreground["symbol"]) & background_genes
)

print("Mapped background genes:", len(background_genes))
print("Mapped significant genes:", len(foreground_genes))

if not background_genes or not foreground_genes:
    raise ValueError("No genes available for enrichment.")

coverage = len(foreground) / len(significant)

print(f"Significant gene mapping coverage: {coverage:.2%}")

if coverage < 0.80:
    raise ValueError(
        "Gene mapping coverage below 80%. "
        "Review mapping before enrichment."
    )


# ============================================================
# 5. PATHWAY ENRICHMENT
# ============================================================

print("\n=== PATHWAY ENRICHMENT ===")

libraries = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human"
]

all_results = []

for library in libraries:

    print(f"\nAnalyzing: {library}")

    gene_sets = gp.get_library(
        name=library,
        organism="Human"
    )

    universe = background_genes
    selected = foreground_genes

    N = len(universe)
    n = len(selected)

    records = []

    for term, members in gene_sets.items():

        members = set(members) & universe

        M = len(members)

        if M == 0:
            continue

        overlap = selected & members
        k = len(overlap)

        pvalue = hypergeom.sf(
            k - 1,
            N,
            M,
            n
        )

        records.append({
            "Term": term,
            "P-value": pvalue,
            "Overlap": k,
            "Term_size": M,
            "Genes": ";".join(sorted(overlap)),
            "Library": library
        })

    result = pd.DataFrame(records)

    if result.empty:
        print("No pathways found.")
        continue

    result["Adjusted P-value"] = multipletests(
        result["P-value"],
        method="fdr_bh"
    )[1]

    result = result.sort_values(
        "Adjusted P-value"
    )

    all_results.append(result)

    output = TABLES / f"{library}_enrichment.csv"

    result.to_csv(
        output,
        index=False
    )

    significant_pathways = result[
        result["Adjusted P-value"] < 0.05
    ]

    print("Terms tested:", len(result))
    print("Significant pathways:", len(significant_pathways))
    print("Saved:", output)


# ============================================================
# 6. COMBINE RESULTS
# ============================================================

print("\n=== SAVING COMBINED RESULTS ===")

if not all_results:
    raise ValueError("No enrichment results generated.")

combined = pd.concat(
    all_results,
    ignore_index=True
)

combined_output = (
    TABLES / "tcga_brca_pathway_enrichment.csv"
)

combined.to_csv(
    combined_output,
    index=False
)

print("Total pathway records:", len(combined))
print("Saved:", combined_output)


# ============================================================
# 7. FINAL SUMMARY
# ============================================================

print("\n================================")
print("PATHWAY ENRICHMENT SUMMARY")
print("================================")

print("Genes tested:", len(df))
print("Significant DEGs:", len(significant))
print("Mapped background genes:", len(background_genes))
print("Mapped significant genes:", len(foreground_genes))

for library in libraries:

    subset = combined[
        combined["Library"] == library
    ]

    count = (
        subset["Adjusted P-value"] < 0.05
    ).sum()

    print(f"{library}: {count} significant pathways")

print("\nSTATUS: PASS")
print("PATHWAY ENRICHMENT COMPLETED")
