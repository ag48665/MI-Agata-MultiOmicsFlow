from pathlib import Path
import hashlib
import shutil
import zipfile

import pandas as pd


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUTS_DIR = (
    PROJECT_ROOT
    / "outputs"
)

SUBMISSION_DIR = (
    OUTPUTS_DIR
    / "submission"
)

PUBLICATION_DIR = (
    OUTPUTS_DIR
    / "publication_package"
)

PACKAGE_CONTENT_DIR = (
    PUBLICATION_DIR
    / "package"
)

ZIP_FILE = (
    PUBLICATION_DIR
    / "MI_Agata_MultiOmics_Month2_Publication_Package.zip"
)

MANIFEST_FILE = (
    PUBLICATION_DIR
    / "publication_manifest.csv"
)

CHECKSUM_FILE = (
    PUBLICATION_DIR
    / "publication_package_sha256.txt"
)

README_FILE = (
    PUBLICATION_DIR
    / "README.md"
)


# =====================================================
# REQUIRED FILES
# =====================================================

REQUIRED_SUBMISSION_FILES = [
    "MI_Agata_MultiOmics_Month2_Final_Manuscript.docx",
    "MI_Agata_MultiOmics_Month2_Final_Manuscript.pdf",
    "figure_captions.md",
    "table_captions.md",
    "manuscript_outline.md",
    "submission_manifest.csv",
    "submission_qc.csv",
    "submission_qc_report.md",
    "manuscript_export_qc.txt",
]


# =====================================================
# HELPERS
# =====================================================

def sha256_file(
    file_path: Path,
) -> str:
    """Calculate SHA256 checksum."""

    digest = hashlib.sha256()

    with file_path.open(
        "rb"
    ) as handle:

        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(
                chunk
            )

    return digest.hexdigest()


def validate_submission() -> None:
    """Check that essential submission files exist."""

    if not SUBMISSION_DIR.exists():
        raise FileNotFoundError(
            f"Submission directory not found:\n"
            f"{SUBMISSION_DIR}"
        )

    missing = []

    for filename in REQUIRED_SUBMISSION_FILES:

        file_path = (
            SUBMISSION_DIR
            / filename
        )

        if (
            not file_path.exists()
            or file_path.stat().st_size <= 0
        ):
            missing.append(
                filename
            )

    if missing:
        raise FileNotFoundError(
            "Missing or empty required submission files:\n"
            + "\n".join(
                missing
            )
        )


def reset_publication_directory() -> None:
    """Rebuild final publication directory."""

    if PUBLICATION_DIR.exists():
        shutil.rmtree(
            PUBLICATION_DIR
        )

    PUBLICATION_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    PACKAGE_CONTENT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


def copy_submission_tree() -> int:
    """Copy complete submission folder into package."""

    file_count = 0

    for source in sorted(
        SUBMISSION_DIR.rglob("*")
    ):

        relative = source.relative_to(
            SUBMISSION_DIR
        )

        destination = (
            PACKAGE_CONTENT_DIR
            / relative
        )

        if source.is_dir():

            destination.mkdir(
                parents=True,
                exist_ok=True,
            )

            continue

        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copy2(
            source,
            destination,
        )

        file_count += 1

    return file_count


# =====================================================
# MANIFEST
# =====================================================

def create_manifest() -> pd.DataFrame:
    """Create manifest of publication package contents."""

    rows = []

    for file_path in sorted(
        PACKAGE_CONTENT_DIR.rglob("*")
    ):

        if not file_path.is_file():
            continue

        relative = (
            file_path
            .relative_to(
                PACKAGE_CONTENT_DIR
            )
            .as_posix()
        )

        rows.append(
            {
                "file":
                    relative,
                "size_bytes":
                    file_path.stat().st_size,
                "sha256":
                    sha256_file(
                        file_path
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =====================================================
# README
# =====================================================

def create_readme(
    file_count: int,
) -> str:
    """Create publication package README."""

    lines = [
        "# MI-Agata-MultiOmicsFlow",
        "",
        "## Month 2 Final Publication Package",
        "",
        (
            "This directory contains the final publication-oriented "
            "outputs generated from the Month 2 clinical-genomic-"
            "epigenomic workflow."
        ),
        "",
        "## Package contents",
        "",
        (
            f"- Total packaged files: {file_count}"
        ),
        "",
        "The package includes:",
        "",
        "- final Word manuscript;",
        "- final PDF manuscript;",
        "- main figures;",
        "- main tables;",
        "- supplementary figures;",
        "- supplementary tables;",
        "- figure captions;",
        "- table captions;",
        "- manuscript outline;",
        "- submission QC results;",
        "- manuscript export QC;",
        "- submission manifest.",
        "",
        "## Final manuscript",
        "",
        (
            "`package/"
            "MI_Agata_MultiOmics_Month2_Final_Manuscript.docx`"
        ),
        "",
        (
            "`package/"
            "MI_Agata_MultiOmics_Month2_Final_Manuscript.pdf`"
        ),
        "",
        "## Integrity",
        "",
        (
            "`publication_manifest.csv` contains SHA256 checksums "
            "for every packaged file."
        ),
        "",
        (
            "`publication_package_sha256.txt` contains the SHA256 "
            "checksum of the final ZIP archive."
        ),
        "",
        "## Validation status",
        "",
        (
            "The analytical Month 2 pipeline passed its automated "
            "final QC audit."
        ),
        "",
        (
            "The submission-material package also passed its "
            "automated QC audit."
        ),
        "",
        "## Important limitation",
        "",
        (
            "The predictive model remains a development-stage "
            "research model based on internal validation."
        ),
        "",
        (
            "Independent external validation is required before "
            "clinical interpretation or deployment."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


# =====================================================
# ZIP
# =====================================================

def create_zip() -> None:
    """Create final publication ZIP archive."""

    print(
        "Creating final ZIP archive..."
    )

    with zipfile.ZipFile(
        ZIP_FILE,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:

        # Add package contents
        for file_path in sorted(
            PACKAGE_CONTENT_DIR.rglob("*")
        ):

            if not file_path.is_file():
                continue

            archive_name = (
                Path("package")
                / file_path.relative_to(
                    PACKAGE_CONTENT_DIR
                )
            )

            archive.write(
                file_path,
                arcname=str(
                    archive_name
                ),
            )

        # Add top-level publication files
        archive.write(
            README_FILE,
            arcname="README.md",
        )

        archive.write(
            MANIFEST_FILE,
            arcname="publication_manifest.csv",
        )


# =====================================================
# ZIP VALIDATION
# =====================================================

def validate_zip() -> tuple[
    int,
    int,
]:
    """Validate ZIP structure and integrity."""

    if not ZIP_FILE.exists():
        raise FileNotFoundError(
            f"ZIP file not created:\n{ZIP_FILE}"
        )

    if ZIP_FILE.stat().st_size <= 0:
        raise ValueError(
            "ZIP file was created but is empty."
        )

    with zipfile.ZipFile(
        ZIP_FILE,
        mode="r",
    ) as archive:

        bad_file = archive.testzip()

        if bad_file is not None:
            raise RuntimeError(
                f"Corrupted ZIP member: {bad_file}"
            )

        members = archive.namelist()

    return (
        len(members),
        ZIP_FILE.stat().st_size,
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Build final publication package."""

    print(
        "Validating submission materials..."
    )

    validate_submission()

    print(
        "Preparing publication package directory..."
    )

    reset_publication_directory()

    print(
        "Copying submission materials..."
    )

    copied_files = (
        copy_submission_tree()
    )

    print(
        f"Copied files: "
        f"{copied_files}"
    )

    # =================================================
    # MANIFEST
    # =================================================

    print(
        "Creating publication manifest..."
    )

    manifest = create_manifest()

    manifest.to_csv(
        MANIFEST_FILE,
        index=False,
    )

    # =================================================
    # README
    # =================================================

    README_FILE.write_text(
        create_readme(
            file_count=(
                len(manifest)
            )
        ),
        encoding="utf-8",
    )

    # =================================================
    # ZIP
    # =================================================

    create_zip()

    # =================================================
    # ZIP QC
    # =================================================

    print(
        "Validating final ZIP archive..."
    )

    member_count, (
        zip_size
    ) = validate_zip()

    zip_sha256 = sha256_file(
        ZIP_FILE
    )

    CHECKSUM_FILE.write_text(
        (
            f"{zip_sha256}  "
            f"{ZIP_FILE.name}\n"
        ),
        encoding="utf-8",
    )

    # =================================================
    # DONE
    # =================================================

    print(
        "\n========================================"
    )

    print(
        "PUBLICATION PACKAGE COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        f"\nFiles in package directory: "
        f"{len(manifest)}"
    )

    print(
        f"Files inside ZIP: "
        f"{member_count}"
    )

    print(
        f"ZIP size: "
        f"{zip_size:,} bytes"
    )

    print(
        "\nZIP SHA256:"
    )

    print(
        zip_sha256
    )

    print(
        "\nPublication directory:\n"
        f"{PUBLICATION_DIR}"
    )

    print(
        "\nFinal ZIP:\n"
        f"{ZIP_FILE}"
    )

    print(
        "\nManifest:\n"
        f"{MANIFEST_FILE}"
    )

    print(
        "\nChecksum file:\n"
        f"{CHECKSUM_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()