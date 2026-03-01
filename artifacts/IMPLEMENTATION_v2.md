# Implementation Plan v2 — Annotations, Legends, HSP90 Explanation

Based on resolved design questions in `ISSUES_v2.md`.

---

## Step 1 — Gene Extractor: Separate annotation loading from gene lookup

**File:** `crossgene/gene_extractor.py`

### Changes:
1. **Add `load_features()` function** — loads sub-gene features from an annotation GTF:
   ```python
   def load_features(
       gene: GeneRecord,
       annotation_gtf_path: str,
       transcript_mode: str = "canonical",  # "canonical" or "all"
       feature_types: set[str] = {"exon", "CDS"},
   ) -> GeneRecord:
   ```
   - Read annotation GTF with gtfparse
   - Filter rows by `gene_id == gene.gene_id`
   - If `transcript_mode == "canonical"`: filter transcript rows for `tag` containing `"Ensembl_canonical"`, then get that transcript_id and filter features by it
   - If `transcript_mode == "all"`: use all features for this gene_id (deduplicating by coordinates)
   - Filter by `feature_type in feature_types`
   - Convert GTF 1-based inclusive → 0-based half-open
   - Return new GeneRecord with updated `features` list

2. **Keep existing `lookup_gene()`** unchanged — it still uses `genes.gtf` for fast gene lookup (gene-level rows only)

### Rationale:
- Separation of concerns: `genes.gtf` (small, fast) for gene lookup; `mainChr.gtf` (large, detailed) for annotation
- The `mainChr.gtf` is ~2GB, so we only read it when annotation is needed (when "plot" format is requested)
- gtfparse caches DataFrames, so reading the annotation GTF once for both genes is efficient

### Checkpoint:
- Load HSP90AA1 features from `mainChr.gtf`, verify exon and CDS features are populated
- Verify canonical transcript filter works (should select HSP90AA1-201 with `Ensembl_canonical` tag)

---

## Step 2 — CLI: Add new parameters and wire annotation loading

**File:** `crossgene/cli.py`

### New CLI parameters:
```
--annotation-gtf    (default: "references/homo_sapiens.109.mainChr.gtf")
--annotation-features  (default: "exon,CDS")
--transcript-mode   (default: "canonical", choices: ["canonical", "all"])
```

### Changes to `main()`:
1. Add the three new Click options
2. Parse `--annotation-features` into a `set[str]`
3. After extracting sequences (`extract_sequence`), if `"plot"` in formats AND `annotation_gtf` exists:
   - Call `load_features()` for both genes with the annotation GTF path, transcript_mode, and feature_types
   - Update `rec_a` and `rec_b` with the returned records
4. Validate `--annotation-gtf` file exists (if provided and plot format requested)

### Checkpoint:
- `crossgene --help` shows the new parameters
- Running with `--output-formats plot` loads features from `mainChr.gtf`

---

## Step 3 — Visualization: Add legend

**File:** `crossgene/visualize.py`

### Changes to `create_circlize_plot()`:
1. After `fig = circos.plotfig()`, add a matplotlib legend at the bottom of the figure
2. Build legend handles dynamically:
   - Always include: same-sense arc (steelblue), antisense arc (firebrick)
   - Conditionally include (only if features present): exon (orange), CDS (forestgreen)
   - Include any other feature types from `--annotation-features` with auto-assigned colors
3. Place legend at bottom center using `fig.legend(..., loc="lower center", ncol=N, bbox_to_anchor=(0.5, -0.02))`
4. Use `matplotlib.patches.Patch` for color swatches

### Color mapping for feature types:
- Extend the existing color constants to a dict:
  ```python
  FEATURE_COLORS = {
      "exon": "orange",
      "CDS": "forestgreen",
      "five_prime_utr": "mediumpurple",
      "three_prime_utr": "goldenrod",
  }
  ```
- Update the feature drawing loop to use this dict (with a fallback color for unknown types)

### Checkpoint:
- Generated PDF has a legend at the bottom
- Legend entries match the features actually drawn (no phantom entries)

---

## Step 4 — Tests

**Files:** `tests/test_gene_extractor.py`, `tests/test_visualize.py`

### New tests for `gene_extractor.py`:
- `test_load_features_canonical`: Verify canonical transcript filtering
- `test_load_features_all`: Verify all-transcripts mode
- `test_load_features_custom_types`: Verify feature_types filter
- `test_load_features_no_annotation`: Verify graceful handling when annotation GTF doesn't exist

### Updated tests for `visualize.py`:
- Verify legend is present in generated figure
- Verify legend entries match features drawn
- Test with features present and without features

### Checkpoint:
- All existing tests still pass
- New tests pass

---

## File Change Summary

| File | Change Type | Description |
|------|------------|-------------|
| `crossgene/gene_extractor.py` | Modified | Add `load_features()` function |
| `crossgene/cli.py` | Modified | Add 3 new CLI params, wire annotation loading |
| `crossgene/visualize.py` | Modified | Add legend, extend feature color mapping |
| `tests/test_gene_extractor.py` | Modified | Add tests for `load_features()` |
| `tests/test_visualize.py` | Modified | Add tests for legend |
| `artifacts/ISSUES_v2.md` | Reference | Design decisions documented |
