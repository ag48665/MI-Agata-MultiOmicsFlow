from pathlib import Path
import shutil
import hashlib
import json
from datetime import datetime


# =====================================
# PATHS
# =====================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUTS_DIR = (
    PROJECT_ROOT
    / "outputs"
)

REPORTS_DIR = (
    OUTPUTS_DIR
    / "reports"
)

TABLES_DIR = (
    OUTPUTS_DIR
    / "tables"
)

FIGURES_DIR = (
    OUTPUTS_DIR
    / "figures"
)

MODELS_DIR = (
    OUTPUTS_DIR
    / "models"
)

PACKAGE_DIR = (
    OUTPUTS_DIR
    / "month2_package"
)

PACKAGE_REPORTS_DIR = (
    PACKAGE_DIR
    / "reports"
)

PACKAGE_TABLES_DIR = (
    PACKAGE_DIR
    / "tables"
)

PACKAGE_FIGURES_DIR = (
    PACKAGE_DIR
    / "figures"
)

PACKAGE_MODELS_DIR = (
    PACKAGE_DIR
    / "models"
)

README_FILE = (
    PACKAGE_DIR
    / "README.md"
)

MANIFEST_FILE = (
    PACKAGE_DIR
    / "manifest.csv"
)

METADATA_FILE = (
    PACKAGE_DIR
    / "package_metadata.json"
)


# =====================================
# FILES TO PACKAGE
# =====================================

REPORT_FILES = [
    "methylation_data_overview.md",
    "methylation_qc_report.md",
    "methylation_patient_matching.md",
    "epigenomic_integration_report.md",
    "epigenomic_figures_report.md",
    "epigenomic_feature_selection.md",
    "epigenomic_baseline_model_report.md",
    "survival_model_report.md",
    "feature_importance_report.md",
    "month2_final_report.md",
]

TABLE_FILES = [
    "epigenomic_integration_summary.csv",
    "methylation_feature_annotation.csv",
    "methylation_feature_variance.csv",
    "selected_methylation_features.csv",
    "epigenomic_model_comparison.csv",
    "clinical_genomic_classification_report.csv",
    "epigenomic_classification_report.csv",
    "clinical_genomic_feature_importance.csv",
    "epigenomic_feature_importance.csv",
    "cox_model_metrics.csv",
    "cox_model_coefficients.csv",
    "cox_patient_risk_scores.csv",
    "cox_selected_features.csv",
    "combined_feature_importance.csv",
]

FIGURE_FILES = [
    "beta_distribution.png",
    "methylation_pca.png",
    "clinical_genomic_confusion_matrix.png",
    "epigenomic_confusion_matrix.png",
    "clinical_genomic_feature_importance.png",
    "epigenomic_feature_importance.png",
    "cox_risk_groups_survival.png",
    "cox_top_hazard_ratios.png",
    "combined_feature_importance.png",
]

MODEL_FILES = [
    "clinical_genomic_random_forest.pkl",
    "clinical_genomic_epigenomic_random_forest.pkl",
    "epigenomic_cox_model.pkl",
]


# =====================================
# FUNCTIONS
# =====================================

def prepare_package_directories() -> None:
    """Create a clean Month 2 package directory."""

    if PACKAGE_DIR.exists():
        shutil.rmtree(
            PACKAGE_DIR
        )

    for directory in [
        PACKAGE_DIR,
        PACKAGE_REPORTS_DIR,
        PACKAGE_TABLES_DIR,
        PACKAGE_FIGURES_DIR,
        PACKAGE_MODELS_DIR,
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def copy_file_if_exists(
    source: Path,
    destination: Path,
) -> bool:
    """Copy file if it exists."""

    if not source.exists():
        print(
            f"Skipping missing file: "
            f"{source.name}"
        )
        return False

    shutil.copy2(
        source,
        destination,
    )

    print(
        f"Copied: {source.name}"
    )

    return True


def calculate_sha256(
    file_path: Path,
) -> str:
    """Calculate SHA256 checksum."""

    sha256 = hashlib.sha256()

    with file_path.open(
        "rb"
    ) as file_handle:

        for block in iter(
            lambda: file_handle.read(
                1024 * 1024
            ),
            b"",
        ):
            sha256.update(
                block
            )

    return sha256.hexdigest()


def collect_manifest_rows() -> list[dict]:
    """Collect metadata for packaged files."""

    rows = []

    for file_path in sorted(
        PACKAGE_DIR.rglob("*")
    ):

        if not file_path.is_file():
            continue

        if file_path in {
            MANIFEST_FILE,
            METADATA_FILE,
        }:
            continue

        relative_path = (
            file_path
            .relative_to(
                PACKAGE_DIR
            )
        )

        rows.append(
            {
                "file": (
                    relative_path
                    .as_posix()
                ),
                "size_bytes": (
                    file_path
                    .stat()
                    .st_size
                ),
                "sha256": (
                    calculate_sha256(
                        file_path
                    )
                ),
            }
        )

    return rows


def write_manifest(
    rows: list[dict],
) -> None:
    """Write CSV manifest without pandas."""

    header = (
        "file,size_bytes,sha256\n"
    )

    lines = [
        header
    ]

    for row in rows:

        lines.append(
            f"{row['file']},"
            f"{row['size_bytes']},"
            f"{row['sha256']}\n"
        )

    MANIFEST_FILE.write_text(
        "".join(lines),
        encoding="utf-8",
    )


def create_readme(
    copied_reports: list[str],
    copied_tables: list[str],
    copied_figures: list[str],
    copied_models: list[str],
) -> str:
    """Create package README."""

    lines = [
        "# Month 2 Results Package",
        "",
        "## Epigenomics Integration and AI-Ready Multi-Omics Pipeline Extension",
        "",
        (
            "This package contains the main results generated "
            "during Month 2 of the multi-omics workflow."
        ),
        "",
        "The module extends the Month 1 clinical-genomics pipeline "
        "with DNA methylation data and includes:",
        "",
        "- methylation quality control",
        "- patient-level matching",
        "- clinical-genomic-epigenomic integration",
        "- methylation feature selection",
        "- epigenomic visualization",
        "- Random Forest classification",
        "- Cox proportional hazards survival modelling",
        "- feature-importance analysis",
        "- final Month 2 reporting",
        "",
        "## Package structure",
        "",
        "```text",
        "month2_package/",
        "├── reports/",
        "├── tables/",
        "├── figures/",
        "├── models/",
        "├── README.md",
        "├── manifest.csv",
        "└── package_metadata.json",
        "```",
        "",
        "## Main report",
        "",
        (
            "The main summary is available at:"
        ),
        "",
        "`reports/month2_final_report.md`",
        "",
        "## Reports",
        "",
    ]

    if copied_reports:
        lines.extend(
            [
                f"- `{name}`"
                for name in copied_reports
            ]
        )
    else:
        lines.append(
            "_No reports were packaged._"
        )

    lines.extend(
        [
            "",
            "## Tables",
            "",
        ]
    )

    if copied_tables:
        lines.extend(
            [
                f"- `{name}`"
                for name in copied_tables
            ]
        )
    else:
        lines.append(
            "_No tables were packaged._"
        )

    lines.extend(
        [
            "",
            "## Figures",
            "",
        ]
    )

    if copied_figures:
        lines.extend(
            [
                f"- `{name}`"
                for name in copied_figures
            ]
        )
    else:
        lines.append(
            "_No figures were packaged._"
        )

    lines.extend(
        [
            "",
            "## Models",
            "",
        ]
    )

    if copied_models:
        lines.extend(
            [
                f"- `{name}`"
                for name in copied_models
            ]
        )
    else:
        lines.append(
            "_No models were packaged._"
        )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            (
                "Large raw and processed molecular data files "
                "are intentionally not duplicated in this package."
            ),
            "",
            (
                "The package focuses on reproducible outputs, "
                "model artefacts, figures, summary tables and reports."
            ),
            "",
            (
                "`manifest.csv` contains file sizes and SHA256 "
                "checksums for integrity verification."
            ),
            "",
        ]
    )

    return "\n".join(
        lines
    )


def create_metadata(
    copied_reports: list[str],
    copied_tables: list[str],
    copied_figures: list[str],
    copied_models: list[str],
) -> dict:
    """Create package metadata."""

    return {
        "package_name": (
            "MI-Agata-MultiOmicsFlow Month 2"
        ),
        "module": (
            "Epigenomics Integration and "
            "AI-Ready Multi-Omics Pipeline Extension"
        ),
        "generated_at": (
            datetime.now()
            .astimezone()
            .isoformat()
        ),
        "reports_count": (
            len(copied_reports)
        ),
        "tables_count": (
            len(copied_tables)
        ),
        "figures_count": (
            len(copied_figures)
        ),
        "models_count": (
            len(copied_models)
        ),
        "raw_data_included": False,
        "processed_large_data_included": False,
    }


def copy_group(
    names: list[str],
    source_directory: Path,
    destination_directory: Path,
) -> list[str]:
    """Copy a group of expected output files."""

    copied = []

    for name in names:

        source = (
            source_directory
            / name
        )

        destination = (
            destination_directory
            / name
        )

        if copy_file_if_exists(
            source,
            destination,
        ):
            copied.append(
                name
            )

    return copied


# =====================================
# MAIN
# =====================================

def main() -> None:
    """Package Month 2 outputs."""

    print(
        "Preparing Month 2 results package..."
    )

    prepare_package_directories()

    print(
        "\nPackaging reports..."
    )

    copied_reports = copy_group(
        names=REPORT_FILES,
        source_directory=REPORTS_DIR,
        destination_directory=(
            PACKAGE_REPORTS_DIR
        ),
    )

    print(
        "\nPackaging tables..."
    )

    copied_tables = copy_group(
        names=TABLE_FILES,
        source_directory=TABLES_DIR,
        destination_directory=(
            PACKAGE_TABLES_DIR
        ),
    )

    print(
        "\nPackaging figures..."
    )

    copied_figures = copy_group(
        names=FIGURE_FILES,
        source_directory=FIGURES_DIR,
        destination_directory=(
            PACKAGE_FIGURES_DIR
        ),
    )

    print(
        "\nPackaging models..."
    )

    copied_models = copy_group(
        names=MODEL_FILES,
        source_directory=MODELS_DIR,
        destination_directory=(
            PACKAGE_MODELS_DIR
        ),
    )

    readme = create_readme(
        copied_reports=copied_reports,
        copied_tables=copied_tables,
        copied_figures=copied_figures,
        copied_models=copied_models,
    )

    README_FILE.write_text(
        readme,
        encoding="utf-8",
    )

    metadata = create_metadata(
        copied_reports=copied_reports,
        copied_tables=copied_tables,
        copied_figures=copied_figures,
        copied_models=copied_models,
    )

    METADATA_FILE.write_text(
        json.dumps(
            metadata,
            indent=2,
        ),
        encoding="utf-8",
    )

    manifest_rows = (
        collect_manifest_rows()
    )

    write_manifest(
        manifest_rows
    )

    print(
        "\nFinished successfully."
    )

    print(
        "\nPackage saved to:\n"
        f"{PACKAGE_DIR}"
    )

    print(
        f"\nReports packaged: "
        f"{len(copied_reports):,}"
    )

    print(
        f"Tables packaged: "
        f"{len(copied_tables):,}"
    )

    print(
        f"Figures packaged: "
        f"{len(copied_figures):,}"
    )

    print(
        f"Models packaged: "
        f"{len(copied_models):,}"
    )

    print(
        "\nREADME:\n"
        f"{README_FILE}"
    )

    print(
        "\nManifest:\n"
        f"{MANIFEST_FILE}"
    )


if __name__ == "__main__":
    main()