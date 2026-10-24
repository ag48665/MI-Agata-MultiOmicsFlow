
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
COUNTS_DIR = RAW_DIR / "star_counts"
PROCESSED_DIR = ROOT / "data" / "processed"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

METADATA_FILE = RAW_DIR / "tcga_brca_selected_226.csv"

metadata = pd.read_csv(METADATA_FILE)

required = {
    "file_id", "patient_id", "sample_type"
}
assert required.issubset(metadata.columns)

assert len(metadata) == 226
assert metadata["file_id"].is_unique
assert metadata["patient_id"].nunique() == 113

sample_counts = metadata.groupby(
    ["patient_id", "sample_type"]
).size().unstack(fill_value=0)

assert (sample_counts["Primary Tumor"] == 1).all()
assert (sample_counts["Solid Tissue Normal"] == 1).all()

files = sorted(COUNTS_DIR.glob("*.tsv"))
assert len(files) == 226

counts_dict = {}
sample_rows = []

for i, file in enumerate(files, start=1):
    match = metadata[
        metadata["file_id"].astype(str).apply(
            lambda fid: fid in file.stem
        )
    ]

    if len(match) != 1:
        raise ValueError(f"Metadata mismatch: {file.name}")

    info = match.iloc[0]

    df = pd.read_csv(
        file,
        sep="\t",
        comment="#",
        usecols=["gene_id", "unstranded"]
    )

    df = df[
        df["gene_id"].str.startswith("ENSG")
    ].copy()

    if len(df) != 60660 or not df["gene_id"].is_unique:
        raise ValueError(f"Invalid gene IDs: {file.name}")

    sample_id = str(info["file_id"])

    counts = pd.to_numeric(
        df["unstranded"], errors="raise"
    )

    if counts.isna().any() or (counts < 0).any():
        raise ValueError(f"Invalid counts: {file.name}")

    if not (counts % 1 == 0).all():
        raise ValueError(f"Noninteger counts: {file.name}")

    counts_dict[sample_id] = pd.Series(
        counts.to_numpy(dtype="int64"),
        index=df["gene_id"]
    )

    sample_rows.append({
        "sample_id": sample_id,
        "patient_id": info["patient_id"],
        "sample_type": info["sample_type"],
        "file_name": file.name
    })

    if i % 25 == 0 or i == len(files):
        print(f"Processed {i}/{len(files)}")

matrix = pd.DataFrame(counts_dict)
samples = pd.DataFrame(sample_rows).set_index("sample_id")

assert matrix.shape == (60660, 226)
assert matrix.index.is_unique
assert not matrix.isna().any().any()
assert matrix.columns.equals(samples.index)

matrix.index.name = "gene_id"

matrix.to_csv(
    PROCESSED_DIR / "tcga_brca_raw_counts.csv"
)

samples.to_csv(
    PROCESSED_DIR / "tcga_brca_sample_metadata.csv"
)

print("\nEXPRESSION MATRIX")
print("Genes:", matrix.shape[0])
print("Samples:", matrix.shape[1])
print("Patients:", samples["patient_id"].nunique())
print("Missing values:", int(matrix.isna().sum().sum()))
print("STATUS: PASS")
