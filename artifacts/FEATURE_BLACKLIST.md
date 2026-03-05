# Feature: `--blacklist` Region Filtering

## Summary

Add a `--blacklist <BED>` CLI option that excludes genomic regions from fragment generation. Fragments overlapping any blacklisted region are skipped entirely, resulting in no alignments (and score 0 in BigWig) for those positions.

Typical use case: mask repetitive elements, low-complexity regions, or other problematic regions using a genome-wide BED file (e.g., ENCODE blacklist, RepeatMasker output).

## Design Decisions

| # | Question | Decision |
|---|----------|----------|
| B1 | Where to apply | **Fragment generation** — skip fragments overlapping blacklisted regions |
| B2 | Single or per-gene | **Single file** — one BED file, regions auto-matched by chromosome |
| B3 | Circular plot | **No special indicator** — blacklisted regions simply produce no arcs |
| B4 | BigWig effect | Score 0 naturally (no fragments = no hits = no score) |
| B5 | Multiple files | **Single file only** — user can merge BED files externally |

## Affected Files

| File | Change |
|------|--------|
| `crossgene/cli.py` | Add `--blacklist` option, parse BED, pass blacklist regions to `_run_direction` |
| `crossgene/fragment.py` | Add `blacklist` parameter to `generate_fragments`, skip overlapping fragments |
| `crossgene/bed_parser.py` | Reuse existing `parse_bed` and `filter_and_clip` (no changes needed) |
| `tests/test_fragment.py` | Add tests for blacklist filtering |
| `tests/test_cli.py` | Add test for `--blacklist` argument parsing |

## Implementation Details

### 1. CLI (`cli.py`)

Add a new click option:

```python
@click.option("--blacklist", default=None, type=click.Path(exists=True),
              help="BED file of regions to exclude from fragment generation.")
```

In `main()`:
- If `--blacklist` is provided, parse it once with `parse_bed(blacklist)`.
- For each direction, use `filter_and_clip` to get blacklist regions for the **query** gene's genomic range.
- Pass the filtered blacklist regions to `_run_direction` and then to `generate_fragments`.

### 2. Fragment Generator (`fragment.py`)

Modify `generate_fragments` signature:

```python
def generate_fragments(
    gene: GeneRecord,
    fragment_size: int,
    step_size: int,
    blacklist: list[BedRegion] | None = None,
) -> Path:
```

For each fragment, before writing it to the FASTA:
1. Convert the fragment's sequence-local position to genomic coordinates using `fragment_index_to_genomic`.
2. Check if the genomic interval overlaps any region in `blacklist`.
3. If overlap: skip the fragment (do not write, do not increment `idx`).
4. If no overlap: write as normal.

Overlap check: a fragment `[frag_start, frag_end)` overlaps a blacklist region `[bl_start, bl_end)` if `frag_start < bl_end and frag_end > bl_start`. Since the blacklist is already filtered/clipped to the gene's range, all regions share the same chromosome.

**Performance note:** For step_size=1 on large genes, the overlap check runs per fragment. With a sorted blacklist, we can use a pointer/index approach (advance through blacklist regions as fragment positions increase) to keep this O(n + m) rather than O(n * m). This is worth implementing since genome-wide blacklists can have many regions.

### 3. Logging

Log the number of blacklisted regions overlapping each gene and the number of fragments skipped:

```
Blacklist: 42 regions overlap BRCA1 (chr17:43044295-43170245)
Direction A->B: skipped 1523/5000 fragments due to blacklist
```

## Data Flow

```
BED file
  |
  v
parse_bed() --> all regions
  |
  v
filter_and_clip(regions, query_gene.chrom, query_gene.start, query_gene.end)
  |
  v
blacklist_regions (for this query gene)
  |
  v
generate_fragments(gene, frag_size, step, blacklist=blacklist_regions)
  |  -- for each fragment: check overlap, skip if blacklisted
  v
fragments FASTA (without blacklisted fragments)
  |
  v
align / parse / score (unchanged -- just fewer fragments)
```

## Edge Cases

- **No blacklist regions overlap gene:** All fragments generated as normal. Log "0 regions overlap".
- **All fragments blacklisted:** Empty FASTA, 0 hits, score array all zeros, BigWig all zeros. Existing "No alignments found" warning covers this.
- **Blacklist file empty or all regions on other chromosomes:** No effect. Log accordingly.
- **Fragment partially overlaps blacklist region:** Fragment is still skipped (any overlap = skip).

## Testing

### `tests/test_fragment.py` — new tests

1. **No blacklist** — existing behavior unchanged (pass `blacklist=None`).
2. **Single blacklisted region** — fragments overlapping it are skipped, others kept. Verify fragment count and indices.
3. **Full overlap** — all fragments blacklisted, empty FASTA.
4. **Partial overlap** — fragment edge touches blacklist boundary, verify skip/keep logic.
5. **Minus-strand gene** — verify genomic coordinate conversion is correct for blacklist check.

### `tests/test_cli.py` — new tests

1. **`--blacklist` with valid BED** — verify it's accepted and passed through.
2. **`--blacklist` with nonexistent file** — click should reject (`exists=True`).
