# Implementation Plan: Flanking Regions Feature

## Overview

Add a `--flanking` parameter that extends gene extraction by N bp upstream and downstream, capturing potential regulatory sites. Flanking regions are visually annotated in the circular plot.

## Design Summary

- `GeneRecord` gains `gene_body_start` / `gene_body_end` to track original gene boundaries
- `start` / `end` expand to include flanks (clamped to chromosome bounds)
- All downstream modules (fragment, align, scores, bigwig) work unchanged — they use `start`/`end`
- Visualization draws flanking as a distinct colored band

---

## Step 1: Update `models.py`

Add two optional fields to `GeneRecord`:

```python
gene_body_start: int = -1   # original gene start before flanking (-1 = not set)
gene_body_end: int = -1     # original gene end before flanking (-1 = not set)
```

When `flanking == 0` or unset, these equal `start`/`end`. A helper property or convention: if `gene_body_start == -1`, treat it as `start`.

---

## Step 2: Update `gene_extractor.py`

### 2a: Modify `extract_sequence(gene, genome_path, flanking=0) -> GeneRecord`

- Store original `gene.start` / `gene.end` as `gene_body_start` / `gene_body_end`
- Compute expanded coords:
  ```
  chrom_len = fa.get_reference_length(gene.chrom)
  flanked_start = max(0, gene.start - flanking)
  flanked_end = min(chrom_len, gene.end + flanking)
  ```
- Fetch the expanded region: `fa.fetch(gene.chrom, flanked_start, flanked_end)`
- If minus strand: reverse-complement as before
- Return `GeneRecord` with `start=flanked_start`, `end=flanked_end`, `gene_body_start`, `gene_body_end` set

### 2b: Adjust `load_features()` (no changes expected)

Features use genomic coords filtered by `gene_id`, so they naturally fall within or outside the gene body. Features in flanking regions from other overlapping genes won't appear (filtered by `gene_id`). No code change needed.

---

## Step 3: Update `cli.py`

- Add Click option:
  ```python
  @click.option("--flanking", default=0, show_default=True,
                help="Flanking region size in bp (upstream + downstream)")
  ```
- Pass `flanking` to both `extract_sequence()` calls
- Log the flanking value if > 0

---

## Step 4: Update `visualize.py`

### 4a: Draw flanking region annotations

In the sector drawing loop, after the outer track axis but before feature annotations:

- If `gene.gene_body_start > gene.start` (left flank exists):
  - Draw a rectangle from local 0 to `gene_body_start - gene.start` with a distinct style (e.g., light lavender with subtle edge)
- If `gene.gene_body_end < gene.end` (right flank exists):
  - Draw a rectangle from `gene_body_end - gene.start` to `gene.end - gene.start`
- For minus-strand genes: the local coordinates still work the same way since `start < end` always in genomic coords

### 4b: Add legend entry

- Add "Flanking region" to the legend with the flanking color (only if flanking regions were actually drawn)
- Add flanking color constant: `COLOR_FLANKING = "lavender"` (or similar subtle color)

---

## Step 5: Propagate `gene_body_start`/`gene_body_end` through GeneRecord copies

Check all places where `GeneRecord` is reconstructed (in `load_features`, `extract_sequence`, `lookup_gene`) — ensure the new fields are carried through. `lookup_gene` sets them to `-1` (not yet flanked). `extract_sequence` sets them. `load_features` preserves them from the input record.

---

## Testing

- Run with `--flanking 0` → behavior identical to current (regression check)
- Run with `--flanking 2000` → verify BigWig extends beyond gene body, plot shows flanking bands
- Edge case: gene near chromosome start/end → flanking clamped correctly
- Minus-strand gene → flanking applied correctly in genomic space
