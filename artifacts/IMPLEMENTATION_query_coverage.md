# Implementation Plan: Query Coverage × Identity Scoring

## Problem

The current BigWig scoring uses `identity = num_matches / alignment_block_length`, which only measures how well the *aligned portion* matches. It does not penalize partial alignments — a 50bp fragment where only 30bp aligns at 95% gets a score of 95, even though 40% of the fragment didn't align. In divergent mode this causes nearly all bases to score 85-100%, producing a flat, uninformative BigWig track.

## Solution

Replace raw identity with **effective identity = identity × query_coverage**, where `query_coverage = aligned_length / fragment_length`. This penalizes partial alignments:
- 50bp fragment, 30bp aligned at 95% → 0.95 × (30/50) = **0.57** (was 0.95)
- 50bp fragment, 48bp aligned at 92% → 0.92 × (48/50) = **0.88** (was 0.92)
- 50bp fragment, 50bp aligned at 100% → 1.00 × (50/50) = **1.00** (unchanged)

## Files to modify

| File | Change |
|------|--------|
| `crossgene/models.py` | Add `query_coverage: float` field to `AlignmentHit` |
| `crossgene/parser.py` | Compute query_coverage from PAF fields (query_alignment_end - query_alignment_start) / query_length |
| `crossgene/parser_blast.py` | Compute query_coverage from BLAST fields: alignment_length / qlen |
| `crossgene/scores.py` | Use `hit.identity * hit.query_coverage` instead of `hit.identity` |
| `crossgene/tsv_writer.py` | Add `query_coverage` column to output |

## Step-by-step

### Step 1: Add `query_coverage` to `AlignmentHit` (models.py)

Add field after `identity`:
```python
query_coverage: float  # 0.0 - 1.0, fraction of query fragment that aligned
```

### Step 2: Compute query_coverage in PAF parser (parser.py)

PAF columns already available:
- Column 1: `query_length` (total length of query sequence)
- Column 2: `query_start` (0-based start of alignment in query)
- Column 3: `query_end` (end of alignment in query)

Compute: `query_coverage = (query_end - query_start) / query_length`

Currently these columns are not parsed — need to extract `query_length` (col 1), `query_start` (col 2), `query_end` (col 3) and compute coverage.

### Step 3: Compute query_coverage in BLAST parser (parser_blast.py)

BLAST tabular format already includes:
- `length` (col index 3): alignment length
- `qlen` (col index 12): query sequence length

Compute: `query_coverage = length / qlen`

Need to add `_COL_QLEN = 12` constant and parse it.

### Step 4: Update scoring (scores.py)

In `compute_scores()`, change:
```python
# Before:
if hit.identity > best_per_fragment[key]:
    best_per_fragment[key] = hit.identity

# After:
effective_identity = hit.identity * hit.query_coverage
if effective_identity > best_per_fragment[key]:
    best_per_fragment[key] = effective_identity
```

### Step 5: Add query_coverage to TSV output (tsv_writer.py)

Add `"query_coverage"` to `COLUMNS` list and include `f"{h.query_coverage:.4f}"` in the row output.

## Impact

- BigWig scores will now range more widely (especially in divergent mode), with partial alignments scored lower
- TSV output gains one additional column
- No CLI changes needed — this changes the internal metric, not the interface
- Existing `--min-quality` filter still applies to raw identity (not effective identity), preserving its documented behavior as an identity threshold
