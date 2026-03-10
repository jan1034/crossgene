# Hit BED Export — Implementation Plan

Based on decisions in [FEATURE_hit_bed_export.md](./FEATURE_hit_bed_export.md).

---

## Steps

### Step 1 — `write_hit_beds()` Function (`bed_writer.py`)

Create a new module that takes a list of `AlignmentHit` and writes two BED9 files: one for the query side and one for the target side, with numbered cross-references.

**File:** `crossgene/bed_writer.py`

**Function:**

```python
def write_hit_beds(
    hits: list[AlignmentHit],
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    direction_label: str,       # "AtoB" or "BtoA"
    output_dir: str,
    primary_only: bool = True,
) -> tuple[str, str]:
    """Write two BED9 files with numbered cross-references.

    Returns (query_bed_path, target_bed_path).
    """
```

**Algorithm:**

1. **Filter:** If `primary_only`, keep only hits where `h.is_primary is True`.
2. **Determine prefixes** from `direction_label`:
   - `"AtoB"` → query prefix `"A"`, target prefix `"B"`
   - `"BtoA"` → query prefix `"B"`, target prefix `"A"`
3. **Assign query-side numbers:** Sort hits by `query_start` ascending (tiebreak: `identity` descending). Enumerate 1..N → build `dict[id(hit), int]` mapping each hit to its query-side number.
4. **Assign target-side numbers:** Sort hits by `target_start` ascending (tiebreak: `identity` descending). Enumerate 1..N → build `dict[id(hit), int]` mapping each hit to its target-side number.
5. **Write query BED:** Iterate hits in query-sorted order. For each hit:
   ```
   chrom    = hit.query_chrom
   start    = hit.query_start
   end      = hit.query_end
   name     = f"{query_prefix}{query_num}-{target_prefix}{target_num}"
   score    = min(int(hit.identity * 1000), 1000)
   strand   = hit.strand
   thickStart = start
   thickEnd   = end
   itemRgb  = "70,130,180" if strand == "+" else "178,34,34"
   ```
6. **Write target BED:** Iterate hits in target-sorted order. For each hit:
   ```
   chrom    = hit.target_chrom
   start    = hit.target_start
   end      = hit.target_end
   name     = f"{target_prefix}{target_num}-{query_prefix}{query_num}"
   score    = min(int(hit.identity * 1000), 1000)
   strand   = hit.strand
   thickStart = start
   thickEnd   = end
   itemRgb  = "70,130,180" if strand == "+" else "178,34,34"
   ```
7. **Track header** for each file:
   ```
   track name="{gene_name}_hits_{direction_label}" description="{query_gene.name} hits → {target_gene.name}" itemRgb=on
   ```
8. **File naming:**
   ```
   {query_gene.name}_hits_{direction_label}.bed
   {target_gene.name}_hits_{direction_label}.bed
   ```
   Example: `BRCA1_hits_AtoB.bed`, `BRCA2_hits_AtoB.bed`

**Constants:**

```python
COLOR_SENSE = "70,130,180"      # steelblue — matches circlize plot
COLOR_ANTISENSE = "178,34,34"   # firebrick — matches circlize plot
```

**Checkpoint:** Import the module, call with synthetic hits, inspect output files manually.

---

### Step 2 — CLI Integration (`cli.py`)

Wire the new writer into the pipeline.

**Changes to `_parse_formats()`:**

```python
valid = {"bigwig", "tsv", "plot", "bed"}  # add "bed"
```

**New CLI option:**

```python
@click.option(
    "--bed-all-hits", is_flag=True, default=False,
    help="Include secondary alignments in hit BED export (default: primary only).",
)
```

**Changes to `_run_direction()`:**

Add `bed_all_hits` parameter. After the TSV block, add:

```python
if "bed" in formats:
    from crossgene.bed_writer import write_hit_beds

    dir_tag = "AtoB" if direction_label == "A→B" else "BtoA"
    q_bed, t_bed = write_hit_beds(
        hits, query_gene, target_gene, dir_tag, outdir,
        primary_only=not bed_all_hits,
    )
    logger.info("Wrote hit BED: %s, %s", q_bed, t_bed)
```

**Changes to `main()` signature:** Add `bed_all_hits` parameter, pass through to `_run_direction()`.

**Checkpoint:** Run `crossgene --gene-a X --gene-b Y --output-formats bed` and verify 4 BED files are created.

---

### Step 3 — Tests (`tests/test_bed_writer.py`)

**Test cases:**

1. **`test_numbering_and_crossrefs`** — Create 3 synthetic `AlignmentHit` objects with known coordinates. Verify:
   - Query BED entries are sorted by `query_start`
   - Target BED entries are sorted by `target_start`
   - Numbers assigned independently per side (by position)
   - Cross-reference names match: if query file has `A1-B3`, target file has `B3-A1`

2. **`test_bed9_format`** — Verify output has exactly 9 tab-separated columns per line (after track header). Verify:
   - `score` is integer in 0-1000 range
   - `thickStart == start`, `thickEnd == end`
   - `itemRgb` is either `70,130,180` or `178,34,34`

3. **`test_track_header`** — First line starts with `track name=`, contains gene name and direction.

4. **`test_primary_only_filter`** — Create hits with `is_primary=True` and `is_primary=False`. Default call excludes secondary. Call with `primary_only=False` includes all.

5. **`test_strand_coloring`** — `+` strand hits get steelblue RGB, `-` strand hits get firebrick RGB.

6. **`test_empty_hits`** — Empty hit list produces BED files with track header only, no crash.

7. **`test_direction_prefixes`** — `"AtoB"` direction uses A/B prefixes; `"BtoA"` uses B/A prefixes.

**Helper:** Create a `_make_hit()` factory function in the test file to build `AlignmentHit` with sensible defaults and only override needed fields.

**Checkpoint:** `pytest tests/test_bed_writer.py -v` all green.

---

### Step 4 — Polish and Edge Cases

- Verify backward compatibility: existing tests pass without `--output-formats bed`.
- Default `--output-formats` remains `"bigwig,tsv,plot"` — BED is opt-in.
- Test with real gene pair to check files load correctly in IGV.
- Verify BED files are sorted (required for IGV indexing with `igvtools`).
  - Query BED is sorted by `query_start` (guaranteed by algorithm).
  - Target BED is sorted by `target_start` (guaranteed by algorithm).
- Handle edge case: hits with identical `query_start` get different numbers (stable sort by start, then identity descending as tiebreaker ensures deterministic ordering).
- Log hit count in BED output (e.g., "Wrote hit BED: BRCA1_hits_AtoB.bed (142 hits)").

**Checkpoint:** Full end-to-end run with `--output-formats bigwig,tsv,plot,bed`, all outputs correct. All tests pass.

---

## Summary

| Step | Module | Type | Dependencies |
|------|--------|------|-------------|
| 1 | `crossgene/bed_writer.py` | New file | None |
| 2 | `crossgene/cli.py` | Edit | Step 1 |
| 3 | `tests/test_bed_writer.py` | New file | Step 1 |
| 4 | Polish | All | Steps 1–3 |
