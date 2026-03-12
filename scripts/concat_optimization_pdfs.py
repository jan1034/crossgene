#!/usr/bin/env python3
"""Concatenate circlize PDFs from an optimization run into a single multi-page PDF.

Each page is annotated with the sensitivity/stringency parameters used.
Requires pypdf and matplotlib.
"""

import argparse
import io
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from pypdf import PdfReader, PdfWriter


def parse_sens_str(dirname: str) -> tuple[int, int] | None:
    """Extract (sensitivity, stringency) from directory name like 'sens3_str2'."""
    m = re.match(r"sens(\d+)_str(\d+)$", dirname)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None


def make_label_pdf(
    text: str,
    page_width_pt: float,
    page_height_pt: float,
    font_size: int = 14,
) -> io.BytesIO:
    """Create a transparent single-page PDF with a label at the top center."""
    w_in = page_width_pt / 72
    h_in = page_height_pt / 72
    fig, ax = plt.subplots(figsize=(w_in, h_in))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5,
        0.975,
        text,
        transform=ax.transAxes,
        fontsize=font_size,
        ha="center",
        va="top",
        fontfamily="sans-serif",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            alpha=0.85,
        ),
    )
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", transparent=True, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    buf.seek(0)

    # The tight bbox changes the page size — re-render at exact page dimensions
    fig2, ax2 = plt.subplots(figsize=(w_in, h_in))
    fig2.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis("off")
    ax2.text(
        0.5,
        0.975,
        text,
        fontsize=font_size,
        ha="center",
        va="top",
        fontfamily="sans-serif",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            alpha=0.85,
        ),
    )
    buf2 = io.BytesIO()
    fig2.savefig(buf2, format="pdf", transparent=True)
    plt.close(fig2)
    buf2.seek(0)
    return buf2


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate optimization circlize PDFs into one summary PDF."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing sens*_str*/ subdirectories",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output PDF path (default: <input_dir>/optimization_summary.pdf)",
    )
    parser.add_argument(
        "--pdf-name",
        default="*.circlize.pdf",
        help="Glob pattern for PDF filename in each subdirectory (default: *.circlize.pdf)",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title prefix on every page",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        default=14,
        help="Annotation font size in pt (default: 14)",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    output_path = args.output or (input_dir / "optimization_summary.pdf")

    # Discover and sort subdirectories
    entries = []
    for d in input_dir.iterdir():
        if not d.is_dir():
            continue
        parsed = parse_sens_str(d.name)
        if parsed is None:
            continue
        pdfs = list(d.glob(args.pdf_name))
        if not pdfs:
            print(f"Warning: no PDF matching '{args.pdf_name}' in {d.name}, skipping")
            continue
        entries.append((parsed, pdfs[0]))

    if not entries:
        print("Error: no sens*_str*/ directories with matching PDFs found", file=sys.stderr)
        sys.exit(1)

    entries.sort(key=lambda x: x[0])

    writer = PdfWriter()
    merged = 0

    for (sens, stri), pdf_path in entries:
        label = f"Sensitivity: {sens}  |  Stringency: {stri}"
        if args.title:
            label = f"{args.title}\n{label}"

        try:
            reader = PdfReader(pdf_path)
        except Exception as e:
            print(f"Warning: could not read {pdf_path}: {e}", file=sys.stderr)
            continue

        for page in reader.pages:
            # Get page dimensions
            box = page.mediabox
            w = float(box.width)
            h = float(box.height)

            # Create label overlay
            label_buf = make_label_pdf(label, w, h, font_size=args.font_size)
            label_page = PdfReader(label_buf).pages[0]

            # Merge label on top of the original page
            page.merge_page(label_page)
            writer.add_page(page)
            merged += 1

    with open(output_path, "wb") as f:
        writer.write(f)

    total = len(entries)
    print(f"Merged {merged} pages from {total} parameter combinations into {output_path}")


if __name__ == "__main__":
    main()
