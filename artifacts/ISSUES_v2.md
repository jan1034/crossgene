# Issues v2 — Feature Annotations, Legends, and HSP90 Identity

## Issue 1: Gene Structure Annotations from `mainChr.gtf`

### Problem
The current pipeline uses `genes.gtf` for gene lookup, which contains **only gene-level rows** (no exons, CDS, UTR). The `gene_extractor.py` code *does* try to extract sub-features (lines 64-77), but finds none because the GTF has no sub-gene features.

A new reference file `homo_sapiens.109.mainChr.gtf` has been added, which contains **full annotations** (gene, transcript, exon, CDS, UTR, etc.).

### Current Behavior
- `gene_extractor.lookup_gene()` uses `genes.gtf` → finds gene row, finds zero sub-features
- `GeneRecord.features` is an empty list
- Visualization draws no exon/CDS annotations on the gene tracks (despite having the code for it in `visualize.py` lines 95-99)

### Proposed Solution
Add a new CLI parameter `--annotation-gtf` (optional) that points to the full `mainChr.gtf`. When provided, use it to load sub-gene features (exon, CDS, UTR) for the gene records. When not provided, attempt to load features from the main `--gtf` as before (which will yield none for `genes.gtf`).

**Design questions for user:**
- Should we default `--annotation-gtf` to `references/homo_sapiens.109.mainChr.gtf`?
- Which feature types to show? Currently coded: exon, CDS, five_prime_utr, three_prime_utr. Should transcript-level features be included?
- Should there be a `--annotation-features` parameter to select which feature types to display (e.g., `--annotation-features exon,CDS`)?
- The `mainChr.gtf` has multiple transcripts per gene. Should we show all transcripts' features (union), or only the canonical transcript?

### Implementation Impact
- `gene_extractor.py`: Add optional annotation GTF path, load features from it
- `cli.py`: Add `--annotation-gtf` parameter, pass to gene extractor
- `visualize.py`: Already handles features when present (but see issue about which features to display per-track)

---

## Issue 2: Missing Legends in Visualization

### Problem
The circlize plot uses four meaningful colors plus alpha transparency for mapq, but includes **no legend**:
- Steelblue arcs = same-sense alignment (+/+ or -/-)
- Firebrick arcs = antisense alignment
- Orange rectangles = exon annotations
- Forest green rectangles = CDS annotations
- Alpha = mapping quality (transparent = low mapq, opaque = high mapq)

Without a legend, the plot is unreadable for anyone not familiar with the color scheme.

### Proposed Solution
Add a matplotlib legend to the figure after `circos.plotfig()`. The legend should include:
- Same-sense alignment (steelblue)
- Antisense alignment (firebrick)
- Exon (orange) — only if features are present
- CDS (forest green) — only if features are present
- Optionally: a note about alpha = mapping quality

### Implementation Impact
- `visualize.py`: Add legend patches after `circos.plotfig()`, before `fig.savefig()`

---

## Issue 3: HSP90AB1/HSP90AA1 — Expected 86% Identity Not Visible

### Diagnosis
The 86% sequence identity between HSP90AB1 and HSP90AA1 reported in literature is at the **protein/coding sequence (CDS) level**, not at the **genomic DNA level**. The key insight:

| Region | Expected Identity | Fraction of Gene |
|--------|------------------|-----------------|
| CDS (exons) | ~86% (protein), ~80-86% (nucleotide) | ~10-20% of gene body |
| Introns | Much lower (<50% or no alignment) | ~80-90% of gene body |
| UTRs | Moderate | Small fraction |

**Why the tool shows low overall identity:**
1. **The tool compares full genomic sequences** (gene body = exons + introns + UTRs), not just CDS/protein
2. Introns diverge rapidly between paralogs — most intronic fragments will not align above `min_quality=30`
3. Intronic regions dominate the gene body (~80-90%), pulling the average mappability score down
4. With default `fragment_size=50`, some fragments at exon-intron boundaries may partly fail to align

**This is not a bug** — it is the expected behavior for a genomic-level comparison tool. The exonic regions *should* show peaks of ~80-86% identity in the BigWig track, while intronic regions show near-zero.

### How Issue 1 Helps
Once gene structure annotations are added (Issue 1), the user will be able to **visually correlate** high-identity peaks with exon positions in:
- The BigWig track in IGV (with exon annotations overlay)
- The circlize plot (exon/CDS colored bars on the gene track, with arcs connecting matching regions)

### Additional Suggestions (for discussion)
- Consider adding a **CDS-only mode** (`--cds-only`) that extracts and compares only coding sequences, which would produce the expected ~86% identity
- Consider adding an **exon-only mode** for mRNA-level comparison
- These modes would concatenate exonic/CDS sequences before fragmentation

---

## Summary of Changes

| Issue | Files Affected | Complexity |
|-------|---------------|------------|
| 1. Annotation GTF | `cli.py`, `gene_extractor.py`, (`visualize.py` already handles it) | Medium |
| 2. Legends | `visualize.py` | Low |
| 3. HSP90 identity | No code change needed — it's expected behavior. But: document it, and issues 1+2 make it understandable visually | N/A (explanation) |

## Resolved Design Questions

| # | Question | Decision |
|---|----------|----------|
| 1 | `--annotation-gtf` default | Default to `references/homo_sapiens.109.mainChr.gtf` |
| 2 | Feature type selection | Add `--annotation-features` param (default: `exon,CDS`) |
| 3 | Multiple transcripts | Add `--transcript-mode` param with values `canonical` (default) / `all`. Canonical = filter by `tag "Ensembl_canonical"` in GTF |
| 4 | CDS-only mode | No — not implementing |
| 5 | Legend placement | Bottom of figure |
