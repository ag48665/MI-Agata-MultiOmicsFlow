from pathlib import Path
import hashlib

import pandas as pd


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

SUBMISSION_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "submission"
)

MAIN_FIGURES_DIR = (
    SUBMISSION_DIR
    / "main_figures"
)

MAIN_TABLES_DIR = (
    SUBMISSION_DIR
    / "main_tables"
)

SUPPLEMENTARY_DIR = (
    SUBMISSION_DIR
    / "supplementary"
)

SUPPLEMENTARY_FIGURES_DIR = (
    SUPPLEMENTARY_DIR
    / "figures"
)

SUPPLEMENTARY_TABLES_DIR = (
    SUPPLEMENTARY_DIR
    / "tables"
)

FIGURE_CAPTIONS_FILE = (
    SUBMISSION_DIR
    / "figure_captions.md"
)

TABLE_CAPTIONS_FILE = (
    SUBMISSION_DIR
    / "table_captions.md"
)

MANUSCRIPT_OUTLINE_FILE = (
    SUBMISSION_DIR
    / "manuscript_outline.md"
)

README_FILE = (
    SUBMISSION_DIR
    / "README.md"
)

MANIFEST_FILE = (
    SUBMISSION_DIR
    / "submission_manifest.csv"
)

QC_TABLE_FILE = (
    SUBMISSION_DIR
    / "submission_qc.csv"
)

QC_REPORT_FILE = (
    SUBMISSION_DIR
    / "submission_qc_report.md"
)


# =====================================================
# EXPECTED FILES
# =====================================================

EXPECTED_MAIN_FIGURES = [
    "Figure_1_model_validation.png",
    "Figure_2_precision_recall.png",
    "Figure_3_calibration.png",
    "Figure_4_confusion_matrix.png",
    "Figure_5_survival.png",
    "Figure_6_feature_importance.png",
]

EXPECTED_MAIN_TABLES = [
    "Table_1_model_validation.csv",
    "Table_2_threshold_optimization.csv",
    "Table_3_calibration.csv",
    "Table_4_final_model_metrics.csv",
    "Table_5_key_results.csv",
]

EXPECTED_SUPPORT_FILES = [
    "figure_captions.md",
    "table_captions.md",
    "manuscript_outline.md",
    "README.md",
    "submission_manifest.csv",
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
    """Append one QC result."""

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


def check_file_exists(
    checks: list[dict],
    file_path: Path,
    category: str,
) -> None:
    """Check file existence and size."""

    if not file_path.exists():

        add_check(
            checks,
            category,
            file_path.name,
            "FAIL",
            f"Missing: {file_path}",
        )

        return

    size = file_path.stat().st_size

    if size <= 0:

        add_check(
            checks,
            category,
            file_path.name,
            "FAIL",
            "File exists but is empty",
        )

        return

    add_check(
        checks,
        category,
        file_path.name,
        "PASS",
        f"Exists; size={size:,} bytes",
    )


# =====================================================
# STRUCTURE QC
# =====================================================

def check_directory_structure(
    checks: list[dict],
) -> None:
    """Check expected submission directories."""

    directories = [
        SUBMISSION_DIR,
        MAIN_FIGURES_DIR,
        MAIN_TABLES_DIR,
        SUPPLEMENTARY_DIR,
        SUPPLEMENTARY_FIGURES_DIR,
        SUPPLEMENTARY_TABLES_DIR,
    ]

    for directory in directories:

        if directory.exists():

            add_check(
                checks,
                "directory",
                directory.name,
                "PASS",
                str(directory),
            )

        else:

            add_check(
                checks,
                "directory",
                directory.name,
                "FAIL",
                f"Missing directory: {directory}",
            )


# =====================================================
# MAIN FILES QC
# =====================================================

def check_main_figures(
    checks: list[dict],
) -> None:
    """Check all expected main figures."""

    for filename in EXPECTED_MAIN_FIGURES:

        check_file_exists(
            checks,
            MAIN_FIGURES_DIR / filename,
            "main_figure",
        )


def check_main_tables(
    checks: list[dict],
) -> None:
    """Check all expected main tables."""

    for filename in EXPECTED_MAIN_TABLES:

        check_file_exists(
            checks,
            MAIN_TABLES_DIR / filename,
            "main_table",
        )


def check_support_files(
    checks: list[dict],
) -> None:
    """Check captions, README, outline and manifest."""

    for filename in EXPECTED_SUPPORT_FILES:

        check_file_exists(
            checks,
            SUBMISSION_DIR / filename,
            "support_file",
        )


# =====================================================
# CSV QC
# =====================================================

def check_csv_tables(
    checks: list[dict],
) -> None:
    """Check that main CSV tables can be read."""

    for filename in EXPECTED_MAIN_TABLES:

        file_path = (
            MAIN_TABLES_DIR
            / filename
        )

        if not file_path.exists():
            continue

        try:

            df = pd.read_csv(
                file_path,
                low_memory=False,
            )

        except Exception as error:

            add_check(
                checks,
                "csv",
                filename,
                "FAIL",
                (
                    f"Could not read CSV: "
                    f"{type(error).__name__}: "
                    f"{error}"
                ),
            )

            continue

        if df.empty:

            add_check(
                checks,
                "csv",
                filename,
                "FAIL",
                "CSV contains no rows",
            )

        else:

            add_check(
                checks,
                "csv",
                filename,
                "PASS",
                (
                    f"Rows={df.shape[0]:,}; "
                    f"columns={df.shape[1]:,}"
                ),
            )


# =====================================================
# CAPTION QC
# =====================================================

def check_figure_captions(
    checks: list[dict],
) -> None:
    """Check that all six figure captions exist."""

    if not FIGURE_CAPTIONS_FILE.exists():
        return

    text = FIGURE_CAPTIONS_FILE.read_text(
        encoding="utf-8"
    )

    for number in range(
        1,
        7,
    ):

        expected = (
            f"Figure {number}"
        )

        if expected in text:

            add_check(
                checks,
                "figure_caption",
                expected,
                "PASS",
                "Caption heading found",
            )

        else:

            add_check(
                checks,
                "figure_caption",
                expected,
                "FAIL",
                "Caption heading missing",
            )


def check_table_captions(
    checks: list[dict],
) -> None:
    """Check that all five table captions exist."""

    if not TABLE_CAPTIONS_FILE.exists():
        return

    text = TABLE_CAPTIONS_FILE.read_text(
        encoding="utf-8"
    )

    for number in range(
        1,
        6,
    ):

        expected = (
            f"Table {number}"
        )

        if expected in text:

            add_check(
                checks,
                "table_caption",
                expected,
                "PASS",
                "Caption heading found",
            )

        else:

            add_check(
                checks,
                "table_caption",
                expected,
                "FAIL",
                "Caption heading missing",
            )


# =====================================================
# MANIFEST QC
# =====================================================

def check_manifest(
    checks: list[dict],
) -> None:
    """Validate manifest and checksums."""

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

    missing_columns = (
        required_columns
        - set(
            manifest.columns
        )
    )

    if missing_columns:

        add_check(
            checks,
            "manifest",
            "Required columns",
            "FAIL",
            (
                f"Missing columns: "
                f"{sorted(missing_columns)}"
            ),
        )

        return

    add_check(
        checks,
        "manifest",
        "Required columns",
        "PASS",
        "file, size_bytes, sha256",
    )

    missing_files = 0
    size_mismatches = 0
    hash_mismatches = 0

    for _, row in (
        manifest.iterrows()
    ):

        relative_path = Path(
            str(
                row["file"]
            )
        )

        file_path = (
            SUBMISSION_DIR
            / relative_path
        )

        if not file_path.exists():

            missing_files += 1
            continue

        actual_size = (
            file_path.stat().st_size
        )

        try:
            expected_size = int(
                row["size_bytes"]
            )
        except (TypeError, ValueError):
            expected_size = -1

        if (
            actual_size
            != expected_size
        ):

            size_mismatches += 1

        actual_hash = sha256_file(
            file_path
        )

        expected_hash = str(
            row["sha256"]
        ).strip()

        if (
            actual_hash
            != expected_hash
        ):

            hash_mismatches += 1

    if missing_files == 0:

        add_check(
            checks,
            "manifest",
            "Manifest files exist",
            "PASS",
            (
                f"{len(manifest):,} "
                "manifest entries checked"
            ),
        )

    else:

        add_check(
            checks,
            "manifest",
            "Manifest files exist",
            "FAIL",
            (
                f"{missing_files} "
                "manifest files missing"
            ),
        )

    if size_mismatches == 0:

        add_check(
            checks,
            "manifest",
            "File sizes",
            "PASS",
            "All sizes match manifest",
        )

    else:

        add_check(
            checks,
            "manifest",
            "File sizes",
            "FAIL",
            (
                f"{size_mismatches} "
                "size mismatches"
            ),
        )

    if hash_mismatches == 0:

        add_check(
            checks,
            "manifest",
            "SHA256 checksums",
            "PASS",
            "All checksums match",
        )

    else:

        add_check(
            checks,
            "manifest",
            "SHA256 checksums",
            "FAIL",
            (
                f"{hash_mismatches} "
                "checksum mismatches"
            ),
        )


# =====================================================
# COUNTS QC
# =====================================================

def check_file_counts(
    checks: list[dict],
) -> None:
    """Check publication material counts."""

    main_figures_count = len(
        list(
            MAIN_FIGURES_DIR.glob(
                "*"
            )
        )
    )

    main_tables_count = len(
        list(
            MAIN_TABLES_DIR.glob(
                "*"
            )
        )
    )

    supplementary_figures_count = len(
        list(
            SUPPLEMENTARY_FIGURES_DIR.glob(
                "*"
            )
        )
    )

    supplementary_tables_count = len(
        list(
            SUPPLEMENTARY_TABLES_DIR.glob(
                "*"
            )
        )
    )

    expected_counts = [
        (
            "Main figures",
            main_figures_count,
            6,
        ),
        (
            "Main tables",
            main_tables_count,
            5,
        ),
        (
            "Supplementary figures",
            supplementary_figures_count,
            6,
        ),
        (
            "Supplementary tables",
            supplementary_tables_count,
            12,
        ),
    ]

    for (
        label,
        actual,
        expected,
    ) in expected_counts:

        if actual == expected:

            add_check(
                checks,
                "count",
                label,
                "PASS",
                f"{actual}",
            )

        else:

            add_check(
                checks,
                "count",
                label,
                "WARNING",
                (
                    f"Expected {expected}; "
                    f"found {actual}"
                ),
            )


# =====================================================
# TEXT QC
# =====================================================

def check_text_files(
    checks: list[dict],
) -> None:
    """Basic text-file QC."""

    text_files = [
        FIGURE_CAPTIONS_FILE,
        TABLE_CAPTIONS_FILE,
        MANUSCRIPT_OUTLINE_FILE,
        README_FILE,
    ]

    for file_path in text_files:

        if not file_path.exists():
            continue

        text = file_path.read_text(
            encoding="utf-8"
        ).strip()

        if len(text) >= 100:

            add_check(
                checks,
                "text",
                file_path.name,
                "PASS",
                (
                    f"{len(text):,} "
                    "characters"
                ),
            )

        else:

            add_check(
                checks,
                "text",
                file_path.name,
                "WARNING",
                (
                    f"Short text file: "
                    f"{len(text):,} characters"
                ),
            )


# =====================================================
# REPORT
# =====================================================

def create_qc_report(
    qc: pd.DataFrame,
    overall_status: str,
) -> str:
    """Create submission QC Markdown report."""

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

    lines = [
        "# Submission Materials QC Report",
        "",
        "## Overall status",
        "",
        f"**{overall_status}**",
        "",
        "## Summary",
        "",
        (
            f"- PASS: "
            f"{counts.get('PASS', 0)}"
        ),
        (
            f"- WARNING: "
            f"{counts.get('WARNING', 0)}"
        ),
        (
            f"- FAIL: "
            f"{counts.get('FAIL', 0)}"
        ),
        "",
    ]

    failures = qc[
        qc["status"]
        == "FAIL"
    ]

    warnings = qc[
        qc["status"]
        == "WARNING"
    ]

    if not failures.empty:

        lines.extend(
            [
                "## Failed checks",
                "",
            ]
        )

        for _, row in (
            failures.iterrows()
        ):

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

        for _, row in (
            warnings.iterrows()
        ):

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
            "## Complete QC table",
            "",
            qc_table,
            "",
            "## Interpretation",
            "",
            (
                "A PASS result confirms that the prepared "
                "submission directory is internally complete "
                "according to the automated checks implemented "
                "in this script."
            ),
            "",
            (
                "This check does not assess journal-specific "
                "formatting requirements or scientific peer "
                "review."
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
    """Run publication-material QC."""

    print(
        "Running submission materials QC..."
    )

    checks = []

    check_directory_structure(
        checks
    )

    check_main_figures(
        checks
    )

    check_main_tables(
        checks
    )

    check_support_files(
        checks
    )

    check_csv_tables(
        checks
    )

    check_figure_captions(
        checks
    )

    check_table_captions(
        checks
    )

    check_manifest(
        checks
    )

    check_file_counts(
        checks
    )

    check_text_files(
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

    report = create_qc_report(
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
        "SUBMISSION MATERIALS QC"
    )

    print(
        "========================================"
    )

    print(
        f"\nSUBMISSION QC STATUS: "
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
        "\nQC table saved to:\n"
        f"{QC_TABLE_FILE}"
    )

    print(
        "\nQC report saved to:\n"
        f"{QC_REPORT_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":
    main()