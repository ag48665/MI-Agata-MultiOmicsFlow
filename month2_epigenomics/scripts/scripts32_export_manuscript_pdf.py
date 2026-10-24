from pathlib import Path
import sys

import win32com.client


# =====================================================
# PATHS
# =====================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

SUBMISSION_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "submission"
)

DOCX_FILE = (
    SUBMISSION_DIR
    / "MI_Agata_MultiOmics_Month2_Final_Manuscript.docx"
)

PDF_FILE = (
    SUBMISSION_DIR
    / "MI_Agata_MultiOmics_Month2_Final_Manuscript.pdf"
)

QC_FILE = (
    SUBMISSION_DIR
    / "manuscript_export_qc.txt"
)


# =====================================================
# WORD CONSTANTS
# =====================================================

WD_EXPORT_FORMAT_PDF = 17

WD_EXPORT_OPTIMIZE_FOR_PRINT = 0

WD_EXPORT_ALL_DOCUMENT = 0

WD_EXPORT_DOCUMENT_CONTENT = 0

WD_EXPORT_CREATE_HEADING_BOOKMARKS = 1


# =====================================================
# HELPERS
# =====================================================

def check_input_file() -> None:
    """Validate source Word manuscript."""

    if not DOCX_FILE.exists():
        raise FileNotFoundError(
            f"Word manuscript not found:\n{DOCX_FILE}"
        )

    if DOCX_FILE.stat().st_size <= 0:
        raise ValueError(
            f"Word manuscript is empty:\n{DOCX_FILE}"
        )


def export_to_pdf() -> None:
    """Export DOCX to PDF using Microsoft Word."""

    print(
        "Starting Microsoft Word..."
    )

    word = None
    document = None

    try:

        word = win32com.client.DispatchEx(
            "Word.Application"
        )

        word.Visible = False

        word.DisplayAlerts = 0

        print(
            "Opening manuscript..."
        )

        document = word.Documents.Open(
            str(
                DOCX_FILE.resolve()
            ),
            ReadOnly=True,
        )

        print(
            "Exporting manuscript to PDF..."
        )

        document.ExportAsFixedFormat(
            OutputFileName=str(
                PDF_FILE.resolve()
            ),
            ExportFormat=(
                WD_EXPORT_FORMAT_PDF
            ),
            OpenAfterExport=False,
            OptimizeFor=(
                WD_EXPORT_OPTIMIZE_FOR_PRINT
            ),
            Range=(
                WD_EXPORT_ALL_DOCUMENT
            ),
            Item=(
                WD_EXPORT_DOCUMENT_CONTENT
            ),
            IncludeDocProps=True,
            KeepIRM=True,
            CreateBookmarks=(
                WD_EXPORT_CREATE_HEADING_BOOKMARKS
            ),
            DocStructureTags=True,
            BitmapMissingFonts=True,
            UseISO19005_1=False,
        )

    finally:

        if document is not None:

            document.Close(
                SaveChanges=False
            )

        if word is not None:

            word.Quit()


def validate_outputs() -> dict:
    """Validate DOCX and PDF after export."""

    results = {
        "docx_exists":
            DOCX_FILE.exists(),

        "pdf_exists":
            PDF_FILE.exists(),

        "docx_size_bytes":
            (
                DOCX_FILE.stat().st_size
                if DOCX_FILE.exists()
                else 0
            ),

        "pdf_size_bytes":
            (
                PDF_FILE.stat().st_size
                if PDF_FILE.exists()
                else 0
            ),
    }

    results[
        "docx_valid"
    ] = (
        results["docx_exists"]
        and results[
            "docx_size_bytes"
        ] > 0
    )

    results[
        "pdf_valid"
    ] = (
        results["pdf_exists"]
        and results[
            "pdf_size_bytes"
        ] > 0
    )

    return results


def create_qc_report(
    results: dict,
) -> str:
    """Create simple export QC report."""

    if (
        results["docx_valid"]
        and results["pdf_valid"]
    ):

        status = "PASS"

    else:

        status = "FAIL"

    lines = [
        "MANUSCRIPT EXPORT QC",
        "=" * 40,
        "",
        f"STATUS: {status}",
        "",
        (
            "DOCX exists: "
            f"{results['docx_exists']}"
        ),
        (
            "DOCX size: "
            f"{results['docx_size_bytes']:,} bytes"
        ),
        "",
        (
            "PDF exists: "
            f"{results['pdf_exists']}"
        ),
        (
            "PDF size: "
            f"{results['pdf_size_bytes']:,} bytes"
        ),
        "",
        f"DOCX path: {DOCX_FILE}",
        "",
        f"PDF path: {PDF_FILE}",
        "",
    ]

    return "\n".join(
        lines
    )


# =====================================================
# MAIN
# =====================================================

def main() -> None:
    """Export final manuscript to PDF."""

    print(
        "Checking final Word manuscript..."
    )

    check_input_file()

    export_to_pdf()

    print(
        "Validating exported files..."
    )

    results = validate_outputs()

    report = create_qc_report(
        results
    )

    QC_FILE.write_text(
        report,
        encoding="utf-8",
    )

    if not results[
        "pdf_valid"
    ]:

        raise RuntimeError(
            "PDF export failed or produced "
            "an empty file."
        )

    print(
        "\n========================================"
    )

    print(
        "FINAL PDF EXPORT COMPLETE"
    )

    print(
        "========================================"
    )

    print(
        "\nDOCX:\n"
        f"{DOCX_FILE}"
    )

    print(
        "\nPDF:\n"
        f"{PDF_FILE}"
    )

    print(
        "\nDOCX size: "
        f"{results['docx_size_bytes']:,} bytes"
    )

    print(
        "PDF size: "
        f"{results['pdf_size_bytes']:,} bytes"
    )

    print(
        "\nExport QC saved to:\n"
        f"{QC_FILE}"
    )

    print(
        "\nFinished successfully."
    )


if __name__ == "__main__":

    try:

        main()

    except Exception as error:

        print(
            "\nPDF EXPORT FAILED"
        )

        print(
            f"{type(error).__name__}: "
            f"{error}"
        )

        sys.exit(1)