from pathlib import Path
import json
import pandas as pd

ROOT = Path("month3_transcriptomics")
BULK = ROOT / "outputs" / "tables"
SC = ROOT / "single_cell" / "outputs" / "tables"
OUT = ROOT / "integration" / "outputs" / "tables"
REPORTS = ROOT / "integration" / "reports"

OUT.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

print("=== BULK + SINGLE-CELL GENE INTEGRATION ===")

mapping = pd.read_csv(BULK / "tcga_brca_gene_mapping.csv")
bulk = pd.read_csv(BULK / "tcga_brca_significant_genes.csv")
markers = pd.read_csv(SC / "leiden_top10_markers.csv")

# Normalize Ensembl IDs by removing version suffixes.
mapping["ensembl_id"] = (
    mapping["ensembl_id"].astype("string")
    .str.split(".").str[0]
)

bulk["ensembl_id"] = (
    bulk["gene_id"].astype("string")
    .str.split(".").str[0]
)

# Remove invalid mapping entries.
mapping["symbol"] = mapping["symbol"].astype("string").str.strip()
mapping = mapping.dropna(subset=["ensembl_id", "symbol"])
mapping = mapping[mapping["symbol"] != ""]
mapping = mapping.drop_duplicates(subset=["ensembl_id", "symbol"])

# Check whether one Ensembl ID maps to multiple symbols.
ambiguous_ids = (
    mapping.groupby("ensembl_id")["symbol"].nunique()
)
ambiguous_ids = ambiguous_ids[ambiguous_ids > 1].index

print("Ambiguous Ensembl IDs:", len(ambiguous_ids))

# Exclude ambiguous mappings rather than choosing arbitrarily.
mapping = mapping[
    ~mapping["ensembl_id"].isin(ambiguous_ids)
].drop_duplicates(subset=["ensembl_id"])

bulk_mapped = bulk.merge(
    mapping[["ensembl_id", "symbol"]],
    on="ensembl_id",
    how="left",
    validate="many_to_one"
)

markers["names"] = markers["names"].astype("string").str.strip()

# Match by exact gene symbol.
overlap = markers.merge(
    bulk_mapped.dropna(subset=["symbol"]),
    left_on="names",
    right_on="symbol",
    how="inner"
)

overlap = overlap.rename(
    columns={
        "group": "leiden_cluster",
        "names": "gene_symbol",
        "scores": "sc_marker_score",
        "logfoldchanges": "sc_log2foldchange",
        "pvals_adj": "sc_padj",
        "log2FoldChange": "bulk_log2foldchange",
        "padj": "bulk_padj",
    }
)

overlap = overlap.sort_values(
    ["leiden_cluster", "sc_marker_score"],
    ascending=[True, False]
)

output_path = OUT / "bulk_single_cell_gene_overlap.csv"
overlap.to_csv(output_path, index=False)

unique_bulk_symbols = set(
    bulk_mapped["symbol"].dropna().astype(str)
)
unique_marker_symbols = set(
    markers["names"].dropna().astype(str)
)
shared_symbols = unique_bulk_symbols & unique_marker_symbols

report = {
    "bulk_significant_rows": int(len(bulk)),
    "bulk_mapped_rows": int(bulk_mapped["symbol"].notna().sum()),
    "bulk_unique_mapped_symbols": int(len(unique_bulk_symbols)),
    "single_cell_marker_rows": int(len(markers)),
    "single_cell_unique_marker_symbols": int(len(unique_marker_symbols)),
    "shared_unique_gene_symbols": int(len(shared_symbols)),
    "overlap_rows": int(len(overlap)),
    "clusters_with_overlap": int(overlap["leiden_cluster"].nunique()),
    "ambiguous_ensembl_ids_excluded": int(len(ambiguous_ids)),
    "note": (
        "Descriptive overlap only. No statistical enrichment "
        "or independent validation is claimed."
    ),
}

report_path = REPORTS / "gene_overlap_report.json"
report_path.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print("\n=== RESULTS ===")
for key, value in report.items():
    print(f"{key}: {value}")

print("\n=== OUTPUT FILES ===")
print(output_path)
print(report_path)
print("\n=== INTEGRATION FINISHED ===")
