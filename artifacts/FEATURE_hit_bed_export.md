# Feature: Hit BED Export with Numbered Cross-References

## Goal

Export alignment hits as BED files so they can be loaded in IGV alongside the BigWig similarity profiles. Each hit is numbered on both the query and target side, with cross-references in the name field so the user can trace which region on Gene A corresponds to which region on Gene B.

## Concept

For each direction (A→B, B→A), produce **two BED files** — one per gene:

```
output/BRCA1_vs_BRCA2/BRCA1_hits_AtoB.bed    # query fragments on Gene A
output/BRCA1_vs_BRCA2/BRCA2_hits_AtoB.bed    # target regions on Gene B
```

### Numbering Scheme

Hits are numbered **independently per gene by genomic position** (sorted ascending). The BED name field encodes a cross-reference to the partner hit number on the other gene:

```
# BRCA1_hits_AtoB.bed  (query side, fragments on Gene A)
chr17  43044295  43044345  A1-B5  92  +
chr17  43044320  43044370  A2-B3  87  +
chr17  43044500  43044550  A3-B1  95  -

# BRCA2_hits_AtoB.bed  (target side, hit regions on Gene B)
chr13  32315480  32315530  B1-A3  95  -
chr13  32890100  32890150  B2-B2  78  +
chr13  32950300  32950350  B3-A2  87  +
chr13  33100000  33100050  B4-B4  65  +
chr13  33200500  33200550  B5-A1  92  +
```

Reading: `A1-B5` means "hit #1 on Gene A maps to hit #5 on Gene B". `B5-A1` is the reciprocal entry on Gene B's BED file.

### Score Field

BED score (column 5) carries the identity percentage (0-100, integer), making hits colorable by score in IGV using the "score" color mode.

### Strand Field

BED strand (column 6) reflects the alignment strand — whether the hit is same-sense (+) or antisense (-) relative to the target.

### Filtering

Only **primary** hits are exported by default to keep the BED files manageable. An optional `--bed-all-hits` flag includes secondary alignments.

### Track Header

Each BED file includes an IGV-compatible track header line:

```
track name="BRCA1_hits_AtoB" description="BRCA1 fragment hits → BRCA2" itemRgb=on
```

### Color by Strand

Use the optional BED9 `itemRgb` column:
- Same-sense (+): steelblue `(70,130,180)` — matches circlize plot colors
- Antisense (-): firebrick `(178,34,34)` — matches circlize plot colors

## Output Format

Full BED9 format:

| Col | Field      | Content                                    |
|-----|------------|--------------------------------------------|
| 1   | chrom      | Chromosome                                 |
| 2   | start      | Hit start (0-based)                        |
| 3   | end        | Hit end (exclusive)                        |
| 4   | name       | `A{n}-B{m}` or `B{m}-A{n}` cross-ref      |
| 5   | score      | Identity * 100 (integer, 0-1000 for BED)   |
| 6   | strand     | Alignment strand (+/-)                     |
| 7   | thickStart | Same as start                              |
| 8   | thickEnd   | Same as end                                |
| 9   | itemRgb    | `70,130,180` (sense) or `178,34,34` (anti) |

Note: BED score is 0-1000 per spec. We map identity (0.0-1.0) to 0-1000.

## Implementation Plan

### Step 1: Add `bed_writer.py` Module

New module: `crossgene/bed_writer.py`

```python
def write_hit_beds(
    hits: list[AlignmentHit],
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    direction_label: str,       # "AtoB" or "BtoA"
    output_dir: str,
    primary_only: bool = True,
) -> tuple[str, str]:
    """Write two BED9 files (query-side, target-side) with numbered cross-refs.

    Returns paths to (query_bed, target_bed).
    """
```

**Algorithm:**

1. Filter hits to primary-only (if `primary_only=True`).
2. Assign **query-side numbers**: sort hits by `query_start` ascending, number 1..N → dict `{hit → query_number}`.
3. Assign **target-side numbers**: sort hits by `target_start` ascending, number 1..M → dict `{hit → target_number}`.
   - Note: N == M (same hits, just sorted differently).
4. Write query BED: iterate hits sorted by `query_start`:
   - `name = f"A{query_number}-B{target_number}"` (or B/A for B→A direction)
   - `score = int(hit.identity * 1000)`
   - `itemRgb` based on strand
5. Write target BED: iterate hits sorted by `target_start`:
   - `name = f"B{target_number}-A{query_number}"` (reciprocal)
   - Same score and strand logic
6. Return both file paths.

**Prefixes:** Use first letter of gene name or just "A"/"B" based on direction to keep names short and unambiguous.

### Step 2: Add `bed` to Output Formats in CLI

In `cli.py`:

1. Add `"bed"` as a recognized output format (alongside `bigwig`, `tsv`, `plot`).
2. Add `--bed-all-hits` flag (default False) — when set, includes secondary alignments.
3. In the per-direction processing loop, after TSV writing, call `write_hit_beds()` if `"bed"` is in output formats.
4. Update the default `--output-formats` to `"bigwig,tsv,plot"` (unchanged — user opts into `bed`).

### Step 3: Tests

In `tests/test_bed_writer.py`:

1. **Test numbering correctness**: Create synthetic hits with known coordinates, verify query-side and target-side numbers are assigned by position order, and cross-references match.
2. **Test BED format**: Verify output is valid BED9 with correct columns, score range 0-1000, track header present.
3. **Test primary filtering**: Verify secondary hits are excluded by default, included with flag.
4. **Test strand coloring**: Verify + hits get steelblue RGB, - hits get firebrick RGB.
5. **Test empty hits**: No crash, produces empty BED files (header only).

## File Changes Summary

| File                      | Change                                      |
|---------------------------|---------------------------------------------|
| `crossgene/bed_writer.py` | **New** — `write_hit_beds()` function       |
| `crossgene/cli.py`        | Add `"bed"` format handling + `--bed-all-hits` |
| `tests/test_bed_writer.py`| **New** — unit tests                        |

## Edge Cases

- **Duplicate coordinates**: Multiple hits with the same query_start get different numbers (stable sort by start, then by identity descending as tiebreaker).
- **B→A direction**: Prefix letters swap — query side uses "B", target side uses "A" — so names remain unambiguous across both direction files.
- **Flanking regions**: Hits in flanking regions are included (they have valid genomic coords). No special handling needed.
- **No hits**: Write BED files with track header only.
