from pathlib import Path
import hashlib
import zipfile

import pandas as pd


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PUBLICATION_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "publication_package"
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

QC_TABLE_FILE = (
    PUBLICATION_DIR
    / "final_release_validation.csv"
)

QC_REPORT_FILE = (
    PUBLICATION_DIR
    / "final_release_validation_report.md"
)


# =====================================================
# EXPECTED ZIP CONTENT
# =====================================================

EXPECTED_FILES = [
    "package/MI_Agata_MultiOmics_Month2_Final_Manuscript.docx",
    "package/MI_Agata_MultiOmics_Month2_Final_Manuscript.pdf",
    "package/figure_captions.md",
    "package/table_captions.md",
    "package/manuscript_outline.md",
    "package/submission_manifest.csv",
    "package/submission_qc.csv",
    "package/submission_qc_report.md",
    "package/manuscript_export_qc.txt",
    "README.md",
    "publication_manifest.csv",
]


# =====================================================
# HELPERS
# =====================================================

def add_check(
    checks: list[dict],
    category: str,
    item: str,
    status: str,
    details: str,
) -> None:
    """Append one release-validation result."""

    checks.append(
        {
            "category": category,
            "item": item,
            "status": status,
            "details": details,
        }
    )


def sha256_file(
    file_path: Path,
) -> str:
    """Calculate SHA256 checksum."""

    digest = hashlib.sha256()

    with file_path.open("rb") as handle:

        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


# =====================================================
# BASIC FILE CHECKS
# =====================================================

def check_release_files(
    checks: list[dict],
) -> None:
    """Check ZIP, manifest and checksum file."""

    for label, file_path in [
        ("ZIP archive", ZIP_FILE),
        ("Publication manifest", MANIFEST_FILE),
        ("ZIP checksum file", CHECKSUM_FILE),
    ]:

        if not file_path.exists():

            add_check(
                checks,
                "release_file",
                label,
                "FAIL",
                f"Missing: {file_path}",
            )

            continue

        size = file_path.stat().st_size

        if size <= 0:

            add_check(
                checks,
                "release_file",
                label,
                "FAIL",
                "File exists but is empty",
            )

        else:

            add_check(
                checks,
                "release_file",
                label,
                "PASS",
                f"Exists; size={size:,} bytes",
            )


# =====================================================
# ZIP INTEGRITY
# =====================================================

def check_zip_integrity(
    checks: list[dict],
) -> list[str]:
    """Validate ZIP and return member names."""

    if not ZIP_FILE.exists():
        return []

    try:

        with zipfile.ZipFile(
            ZIP_FILE,
            mode="r",
        ) as archive:

            bad_file = archive.testzip()

            members = archive.namelist()

            if bad_file is None:

                add_check(
                    checks,
                    "zip",
                    "Archive integrity",
                    "PASS",
                    "No corrupted ZIP members detected",
                )

            else:

                add_check(
                    checks,
                    "zip",
                    "Archive integrity",
                    "FAIL",
                    f"Corrupted member: {bad_file}",
                )

            if len(members) > 0:

                add_check(
                    checks,
                    "zip",
                    "Archive member count",
                    "PASS",
                    f"{len(members):,} members",
                )

            else:

                add_check(
                    checks,
                    "zip",
                    "Archive member count",
                    "FAIL",
                    "ZIP archive contains no files",
                )

            return members

    except zipfile.BadZipFile as error:

        add_check(
            checks,
            "zip",
            "ZIP readable",
            "FAIL",
            f"BadZipFile: {error}",
        )

        return []


# =====================================================
# REQUIRED ZIP CONTENT
# =====================================================

def check_expected_files(
    checks: list[dict],
    members: list[str],
) -> None:
    """Check essential files inside ZIP."""

    member_set = set(members)

    for expected in EXPECTED_FILES:

        if expected in member_set:

            add_check(
                checks,
                "zip_content",
                expected,
                "PASS",
                "Present in ZIP",
            )

        else:

            add_check(
                checks,
                "zip_content",
                expected,
                "FAIL",
                "Missing from ZIP",
            )


# =====================================================
# EMPTY MEMBERS
# =====================================================

def check_zip_member_sizes(
    checks: list[dict],
) -> None:
    """Check that packaged files are not empty."""

    if not ZIP_FILE.exists():
        return

    with zipfile.ZipFile(
        ZIP_FILE,
        mode="r",
    ) as archive:

        empty_files = []

        checked_files = 0

        for info in archive.infolist():

            if info.is_dir():
                continue

            checked_files += 1

            if info.file_size <= 0:
                empty_files.append(
                    info.filename
                )

        if not empty_files:

            add_check(
                checks,
                "zip_content",
                "Non-empty packaged files",
                "PASS",
                f"{checked_files:,} files checked",
            )

        else:

            add_check(
                checks,
                "zip_content",
                "Non-empty packaged files",
                "FAIL",
                (
                    f"{len(empty_files)} empty files: "
                    + ", ".join(empty_files[:10])
                ),
            )


# =====================================================
# MANIFEST CHECK
# =====================================================

def check_manifest(
    checks: list[dict],
) -> None:
    """Validate publication manifest structure."""

    if not MANIFEST_FILE.exists():
        return

    try:

        manifest = pd.read_csv(
            MANIFEST_FILE,
            low_memory=False,
        )

    except Exception as error:

        add_check(
            checks,
            "manifest",
            "Manifest readable",
            "FAIL",
            (
                f"{type(error).__name__}: "
                f"{error}"
            ),
        )

        return

    required_columns = {
        "file",
        "size_bytes",
        "sha256",
    }

    missing = (
        required_columns
        - set(manifest.columns)
    )

    if missing:

        add_check(
            checks,
            "manifest",
            "Manifest columns",
            "FAIL",
            f"Missing columns: {sorted(missing)}",
        )

        return

    add_check(
        checks,
        "manifest",
        "Manifest columns",
        "PASS",
        "file, size_bytes, sha256",
    )

    if manifest.empty:

        add_check(
            checks,
            "manifest",
            "Manifest rows",
            "FAIL",
            "Manifest contains no rows",
        )

    else:

        add_check(
            checks,
            "manifest",
            "Manifest rows",
            "PASS",
            f"{len(manifest):,} entries",
        )

    duplicates = int(
        manifest["file"]
        .duplicated()
        .sum()
    )

    if duplicates == 0:

        add_check(
            checks,
            "manifest",
            "Unique manifest paths",
            "PASS",
            "No duplicated file paths",
        )

    else:

        add_check(
            checks,
            "manifest",
            "Unique manifest paths",
            "FAIL",
            f"{duplicates} duplicated file paths",
        )


# =====================================================
# ZIP SHA256 CHECK
# =====================================================

def check_zip_checksum(
    checks: list[dict],
) -> None:
    """Compare calculated ZIP SHA256 with saved checksum."""

    if (
        not ZIP_FILE.exists()
        or not CHECKSUM_FILE.exists()
    ):
        return

    calculated = sha256_file(
        ZIP_FILE
    )

    checksum_text = (
        CHECKSUM_FILE
        .read_text(
            encoding="utf-8"
        )
        .strip()
    )

    expected = (
        checksum_text
        .split()[0]
        if checksum_text
        else ""
    )

    if not expected:

        add_check(
            checks,
            "checksum",
            "Saved SHA256",
            "FAIL",
            "Checksum file contains no SHA256 value",
        )

        return

    if calculated == expected:

        add_check(
            checks,
            "checksum",
            "ZIP SHA256",
            "PASS",
            calculated,
        )

    else:

        add_check(
            checks,
            "checksum",
            "ZIP SHA256",
            "FAIL",
            (
                f"Expected={expected}; "
                f"calculated={calculated}"
            ),
        )


# =====================================================
# MANUSCRIPT FILES INSIDE ZIP
# =====================================================

def check_manuscript_files(
    checks: list[dict],
) -> None:
    """Check DOCX/PDF sizes inside ZIP."""

    if not ZIP_FILE.exists():
        return

    expected = {
        "package/MI_Agata_MultiOmics_Month2_Final_Manuscript.docx":
            "DOCX manuscript",
        "package/MI_Agata_MultiOmics_Month2_Final_Manuscript.pdf":
            "PDF manuscript",
    }

    with zipfile.ZipFile(
        ZIP_FILE,
        mode="r",
    ) as archive:

        names = set(
            archive.namelist()
        )

        for filename, label in (
            expected.items()
        ):

            if filename not in names:
                continue

            info = archive.getinfo(
                filename
            )

            if info.file_size > 10_000:

                add_check(
                    checks,
                    "manuscript",
                    label,
                    "PASS",
                    (
                        f"Size={info.file_size:,} bytes"
                    ),
                )

            else:

                add_check(
                    checks,
                    "manuscript",
                    label,
                    "WARNING",
                    (
                        f"Unexpectedly small: "
                        f"{info.file_size:,} bytes"
                    ),
                )


# =====================================================
# REPORT
# =====================================================

def create_report(
    qc: pd.DataFrame,
    overall_status: str,
) -> str:
    """Create Markdown final release report."""

    counts = (
        qc["status"]
        .value_counts()
        .to_dict()
    )

    try:

        qc_table = qc.to_markdown(
            index=False
        )

    except ImportError:

        qc_table = (
            "```text\n"
            + qc.to_string(
                index=False
            )
            + "\n```"
        )

    failures = qc[
        qc["status"] == "FAIL"
    ]

    warnings = qc[
        qc["status"] == "WARNING"
    ]

    lines = [
        "# Final Publication Release Validation",
        "",
        "## Overall status",
        "",
        f"**{overall_status}**",
        "",
        "## Summary",
        "",
        f"- PASS: {counts.get('PASS', 0)}",
        f"- WARNING: {counts.get('WARNING', 0)}",
        f"- FAIL: {counts.get('FAIL', 0)}",
        "",
    ]

    if overall_status == "PASS":

        lines.extend(
            [
                "## Release decision",
                "",
                (
                    "The final Month 2 publication package "
                    "passed all automated release-validation checks."
                ),
                "",
                (
                    "The ZIP archive is readable, required "
                    "submission artefacts are present, and the "
                    "saved SHA256 checksum matches the archive."
                ),
                "",
            ]
        )

    elif overall_status == "WARNING":

        lines.extend(
            [
                "## Release decision",
                "",
                (
                    "No critical release failure was detected, "
                    "but warnings should be reviewed before "
                    "external distribution."
                ),
                "",
            ]
        )

    else:

        lines.extend(
            [
                "## Release decision",
                "",
                (
                    "The publication package should not be "
                    "distributed until all failed checks have "
                    "been resolved."
                ),
                "",
            ]
        )

    if not failures.empty:

        lines.extend(
            [
                "## Failed checks",
                "",
            ]
        )

        for _, row in failures.iterrows():

            lines.append(
                (
                    f"- **{row['category']} / "
                    f"{row['item']}**: "
                    f"{row['details']}"
                )
            )

        lines.append("")

    if not warnings.empty:

        lines.extend(
            [
                "## Warnings",
                "",
            ]
        )

        for _, row in warnings.iterrows():

            lines.append(
                (
                    f"- **{row['category']} / "
                    f"{row['item']}**: "
                    f"{row['details']}"
                )
            )

        lines.append("")

    lines.extend(
        [
            "## Complete validation table",
            "",
            qc_table,
            "",
            "## Scientific limitation",
            "",
            (
                "This release validation confirms technical "
                "integrity of the packaged outputs. It does not "
                "constitute independent scientific or clinical "
                "validation of the predictive model."
            ),
            "",
        ]
    )

    return "\n".join(
        lines
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Validate final publication release."""

    print(
        "Running final publication release validation..."
    )

    checks = []

    check_release_files(
        checks
    )

    members = check_zip_integrity(
        checks
    )

    check_expected_files(
        checks,
        members,
    )

    check_zip_member_sizes(
        checks
    )

    check_manifest(
        checks
    )

    check_zip_checksum(
        checks
    )

    check_manuscript_files(
        checks
    )

    qc = pd.DataFrame(
        checks
    )

    if (
        qc["status"] == "FAIL"
    ).any():

        overall_status = "FAIL"

    elif (
        qc["status"] == "WARNING"
    ).any():

        overall_status = "WARNING"

    else:

        overall_status = "PASS"

    qc.to_csv(
        QC_TABLE_FILE,
        index=False,
    )

    report = create_report(
        qc=qc,
        overall_status=overall_status,
    )

    QC_REPORT_FILE.write_text(
        report,
        encoding="utf-8",
    )

    print(
        "\n========================================"
    )

    print(
        "FINAL RELEASE VALIDATION"
    )

    print(
        "========================================"
    )

    print(
        f"\nFINAL RELEASE STATUS: "
        f"{overall_status}"
    )

    print(
        "\nStatus counts:"
    )

    print(
        qc["status"]
        .value_counts()
        .to_string()
    )

    failures = qc[
        qc["status"] == "FAIL"
    ]

    warnings = qc[
        qc["status"] == "WARNING"
    ]

    if not failures.empty:

        print(
            "\nFAILED CHECKS:"
        )

        for _, row in (
            failures.iterrows()
        ):

            print(
                f"- {row['category']} / "
                f"{row['item']}: "
                f"{row['details']}"
            )

    if not warnings.empty:

        print(
            "\nWARNINGS:"
        )

        for _, row in (
            warnings.iterrows()
        ):

            print(
                f"- {row['category']} / "
                f"{row['item']}: "
                f"{row['details']}"
            )

    print(
        "\nValidation table saved to:\n"
        f"{QC_TABLE_FILE}"
    )

    print(
        "\nValidation report saved to:\n"
        f"{QC_REPORT_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()