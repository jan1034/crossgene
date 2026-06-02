# Plan: PDF Concatenation Script for Parameter Optimization

## Problem

After running `process_2x.sh`, each sensitivity x stringency combination produces a separate circlize PDF in its own subdirectory (e.g., `output/actg1-actg2/sens1_str1/Actg1_vs_Actg2.circlize.pdf`). Comparing 25 plots across directories is tedious. Need a single scrollable PDF per optimization run with clear parameter annotations on each page.

## Scope

- Concatenate PDFs **within** one optimization run (one `output/<run>/` directory)
- NOT across run types (e.g., native vs blacklisted stay separate)
- Each page gets a text annotation showing the sensitivity and stringency used
- Output: one combined PDF per run directory

## Existing Structure

```
output/
├── actg1-actg2/                    # native run
│   ├── sens1_str1/
│   │   └── Actg1_vs_Actg2.circlize.pdf
│   ├── sens1_str2/
│   │   └── ...
│   └── sens5_str5/
└── actg1-actg2-no-sines/           # blacklisted run
    ├── sens1_str1/
    └── ...
```

Directory naming convention: `sens{1-5}_str{1-5}/`

## Implementation

### Script: `scripts/concat_optimization_pdfs.py`

Python script using only matplotlib (already in the environment, no new dependencies).

**Approach:** matplotlib's `PdfPages` to create a new multi-page PDF. For each input PDF, render it as an image (via matplotlib's PDF backend reading), then draw it onto a new page with a header annotation.

Since matplotlib can't natively read/embed existing PDFs as images, the practical approach is:

1. Use `matplotlib.backends.backend_pdf.PdfPages` for output
2. Convert each input PDF page to a raster image via matplotlib's `subprocess` call to `pdftoppm` or — more portably — re-render with matplotlib's own PDF reader

**Revised approach — pure matplotlib re-rendering won't work for arbitrary PDFs. Two options:**

### Option A: Use PyPDF2/pypdf (preferred, small dependency)

Install `pypdf` (pure Python, no C deps). Use it to read input PDFs, stamp a text annotation on each page, and merge into one output.

```
pip install pypdf   # ~200 KB, pure Python
```

**Steps:**
1. Glob for `sens*_str*/` directories under the input path
2. Sort by (sensitivity asc, stringency asc) for logical page order
3. Parse `sens{N}_str{M}` from directory name
4. For each: read the circlize PDF, stamp annotation text, append to output
5. Write combined PDF

**Annotation:** Add a text overlay at the top of each page:
- Text: `Sensitivity: {N}  |  Stringency: {M}`
- Semi-transparent white background rectangle so it doesn't clash with the plot
- Font: Helvetica, ~14pt

### Option B: Use matplotlib only (no new deps, but rasterizes)

Convert each PDF to PNG via matplotlib, then compose pages with annotation text. Downside: rasterization loses vector quality.

### Decision: Option A (pypdf)

Vector quality preserved, tiny dependency, clean implementation.

### CLI Interface

```
python scripts/concat_optimization_pdfs.py <input_dir> [options]
```

**Arguments:**
- `input_dir` (positional, required) — path containing `sens*_str*/` subdirectories (e.g., `output/actg1-actg2/`)

**Options:**
- `--output` / `-o` — output PDF path (default: `<input_dir>/optimization_summary.pdf`)
- `--pdf-name` — glob pattern for the PDF filename within each subdirectory (default: `*.circlize.pdf`)
- `--title` — optional title text for the header on every page (e.g., "Actg1 vs Actg2 — native")
- `--font-size` — annotation font size in pt (default: 14)

### Page ordering

Sorted by sensitivity first, then stringency: `(1,1), (1,2), ..., (1,5), (2,1), ..., (5,5)`. This groups all stringency levels for a given sensitivity together, making it easy to see how filtering affects a given alignment.

### Annotation placement

- Top-center of the page
- White background box with slight padding
- Text: `Sensitivity: 1  |  Stringency: 3` (or with optional title prefix)
- Should not overlap the circular plot (the plots have margins at top)

### Implementation Steps

1. Create `scripts/concat_optimization_pdfs.py`
2. Parse args with `argparse` (no click dependency needed for a standalone script)
3. Glob `input_dir/sens*_str*/` and sort
4. For each directory:
   - Extract sens/str from dirname via regex
   - Find the PDF file matching `--pdf-name`
   - Skip with warning if no PDF found
5. Use `pypdf.PdfReader` to read each input PDF
6. Use `pypdf.PdfWriter` for output
7. For each page, use `pypdf` annotations or overlay a small PDF (generated with reportlab or matplotlib) containing the label text
8. Write final merged PDF

### Text overlay without reportlab

`pypdf` alone can't easily add arbitrary text. Two sub-approaches:

**A) Generate a tiny single-page label PDF with matplotlib, then overlay via pypdf's `merge_page`:**
```python
# For each (sens, str) combo:
fig, ax = plt.subplots(figsize=(page_width_in, page_height_in))
ax.text(0.5, 0.97, label_text, transform=ax.transAxes, ...)
ax.axis('off')
fig.savefig(label_buf, format='pdf', transparent=True)

# Then merge:
label_page = PdfReader(label_buf).pages[0]
input_page.merge_page(label_page)
writer.add_page(input_page)
```

This keeps the original plot as vectors and overlays a transparent label PDF on top. matplotlib is already available — no new dependencies beyond `pypdf`.

**B) Full matplotlib rasterize approach (fallback if pypdf not available):**
Read each PDF as an image, draw onto matplotlib figure with annotation. Loses vector quality but needs zero new deps.

### Final design: support both modes

- If `pypdf` is installed: vector merge + matplotlib label overlay (best quality)
- Else: matplotlib-only rasterized fallback with warning about quality loss

### Error handling

- Missing directories: skip with warning
- Missing PDFs in a directory: skip with warning, note in summary
- No PDFs found at all: exit with error
- Print summary: "Merged N/25 pages into output.pdf"

### Testing

Manual — run on existing optimization output and verify:
```bash
conda run -n compare_genes python scripts/concat_optimization_pdfs.py \
    comparisons/mus_musculus/output/actg1-actg2/ \
    --title "Actg1 vs Actg2 (native)"

conda run -n compare_genes python scripts/concat_optimization_pdfs.py \
    comparisons/mus_musculus/output/actg1-actg2-no-sines/ \
    --title "Actg1 vs Actg2 (SINEs blacklisted)"
```
