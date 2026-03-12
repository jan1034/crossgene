# Plan: Remove Reverse-Complement from Sequence Extraction

## Rationale

The current pipeline reverse-complements minus-strand gene sequences in `extract_sequence()`. This complicates every downstream step that needs to map between sequence positions and genomic coordinates (blacklist masking, coordinate translation, etc.) — and already caused a bug in `mask_sequence`. BLAST handles strandedness natively, so the revcomp adds no value.

## What Changes

### 1. `gene_extractor.py` — Remove revcomp

**Lines 220-221:** Delete the minus-strand branch:
```python
# DELETE:
if gene.strand == "-":
    seq = _reverse_complement(seq)
```

The `_reverse_complement()` function can also be removed (lines 15-18) since nothing else uses it.

### 2. `parser_blast.py` — Simplify `_local_to_genomic()`

Currently has a minus-strand branch that reverses the mapping. With genomic-orientation sequences, it becomes the same for both strands:

```python
# BEFORE (strand-dependent):
if gene.strand == "-":
    genomic_end = gene.end - local_start
    genomic_start = gene.end - local_end
else:
    genomic_start = gene.start + local_start
    genomic_end = gene.start + local_end

# AFTER (strand-independent):
genomic_start = gene.start + local_start
genomic_end = gene.start + local_end
```

### 3. `blastn.py` — Simplify `mask_sequence()`

Remove `gene_end` and `strand` parameters. Always use `region.start - gene_start`:

```python
# BEFORE:
def mask_sequence(sequence, regions, gene_start, gene_end, strand):
    ...
    if strand == "-":
        local_start = gene_end - region.end
        local_end = gene_end - region.start
    else:
        local_start = region.start - gene_start
        local_end = region.end - gene_start

# AFTER:
def mask_sequence(sequence, regions, gene_start):
    ...
    local_start = region.start - gene_start
    local_end = region.end - gene_start
```

Update call site in `align_genes()` accordingly.

### 4. `models.py` — Update docstring

Change `sequence` field comment from "sense sequence (revcomp applied if minus strand)" to "genomic-orientation sequence".

### 5. Hit strand semantics — Derive relative strand from both genes

The `strand` field in `AlignmentHit` comes from BLAST's `sstrand` column, which reports whether the *target* aligned on plus or minus strand. Currently, with both sequences revcomp'd to sense orientation, BLAST `sstrand == "plus"` always means same-sense. After removing revcomp, both sequences are in genomic orientation, so `sstrand == "plus"` means "same genomic direction" — which is only same-sense if both genes are on the same strand.

**Resolution:** After removing revcomp, derive the sense/antisense label from both genes' strands:
```python
# In parser_blast.py when building AlignmentHit:
blast_strand = "+" if sstrand == "plus" else "-"
# Both sequences are genomic orientation. BLAST "plus" = same genomic direction.
# Same-sense means same transcription direction, which requires considering both genes.
if query_gene.strand == target_gene.strand:
    relative_strand = blast_strand
else:
    relative_strand = "-" if blast_strand == "+" else "+"
```

Truth table:

| Query strand | Target strand | BLAST sstrand | Same genomic dir? | Same-sense? | Result |
|---|---|---|---|---|---|
| + | + | plus  | yes | yes | + |
| + | + | minus | no  | no  | - |
| - | + | plus  | yes | no  | - |
| - | + | minus | no  | yes | + |
| + | - | plus  | yes | no  | - |
| + | - | minus | no  | yes | + |
| - | - | plus  | yes | yes | + |
| - | - | minus | no  | no  | - |

This keeps the downstream meaning of `hit.strand` unchanged (+ = same-sense, - = antisense).

## Files NOT affected

- **scores.py** — operates on genomic coordinates from `AlignmentHit`, no change
- **bigwig.py** — writes scores at genomic coordinates, no change
- **tsv_writer.py** — passes through `hit.strand`, no change (semantics preserved by fix in #5)
- **bed_writer.py** — uses `hit.strand` for coloring, no change (semantics preserved by fix in #5)
- **visualize.py** — uses `hit.strand` for arc coloring, no change (semantics preserved by fix in #5)

## Cleanup

- Delete `_reverse_complement()` function from `gene_extractor.py` (dead code after removing the revcomp call)
- Remove `gene_end` and `strand` parameters from `mask_sequence()` signature and its call site
- Remove the minus-strand branch from `_local_to_genomic()`
- Update `PLAN_debug_and_sensitivity_refactor.md` to reflect the simplified blacklist flow

## Tests to update

### `tests/test_gene_extractor.py`
- **`TestExtractSequence.test_minus_strand_revcomp`** (line 83): Currently expects `"TTTT" * 25` (revcomp of all A's). After change, should expect `"AAAA" * 25` (raw genomic sequence).
- **`TestReverseComplement`** (line 185): Entire test class tests the deleted `_reverse_complement()` function. **Delete this class.**
- **Line 10 import**: Remove `_reverse_complement` from the import.

### `tests/test_parser_blast.py`
- **`TestCoordinateTranslation.test_query_coords`** (line 86): The minus-strand parametrize case `("-", "51", "250", 1750, 1950)` currently expects reversed genomic coords. After change, should expect `1050, 1250` (same as plus strand — simple `gene_start + local`).
- **`TestCoordinateTranslation.test_target_coords`** (line 100): The minus-strand target case `("-", "101", "300", "plus", 5700, 5900, "+")` currently expects reversed coords. After change, should expect `5100, 5300` (same as plus strand) **and strand `"-"`** because query(+) target(-) with blast "plus" = antisense → `("-", "101", "300", "plus", 5100, 5300, "-")`.

### `tests/test_blastn.py`
- **`TestWriteFastaAndMask.test_mask_sequence`** (line 87): Tests already use the simplified 3-param signature `mask_sequence(seq, regions, gene_start=1000)` — they are currently broken against the old 5-param code and will pass once `mask_sequence()` is simplified. **No change needed.**

### Tests NOT affected
- `test_scores.py`, `test_bigwig.py`, `test_tsv_writer.py`, `test_bed_writer.py`, `test_visualize.py`, `test_bed_parser.py`, `test_cli.py` — no dependency on revcomp behavior.

## Summary of changes

| File | Change | Simplifies? |
|------|--------|-------------|
| gene_extractor.py | Remove revcomp + helper function | Yes |
| parser_blast.py | Remove minus-strand branch in `_local_to_genomic()`, add relative strand derivation | Yes (net) |
| blastn.py | Remove strand/gene_end params from `mask_sequence()` | Yes |
| models.py | Update docstring | Trivial |
| tests | Update expected values | — |
