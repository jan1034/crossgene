# Plan: Direct Gene-vs-Gene BLASTN (Simplified Pipeline)

## Problem Statement

The current fragment→merge pipeline doesn't reliably represent true homologous regions. Fragment-level alignment (50 bp windows) is inherently noisy: short queries produce approximate boundaries, overlapping fragment scores are averaged rather than properly computed, and the merge heuristic chains fragments by adjacency rather than biological signal. Additionally, maintaining two aligners (minimap2 + BLASTN) adds complexity without clear benefit.

## Decision: Direct Gene-vs-Gene BLAST

After discussion, the chosen approach is the simplest possible: **run BLASTN directly with gene A as query against gene B as subject**. BLASTN natively finds all homologous segments (HSPs) with proper scoring, boundaries, and statistics. No fragmentation, no merging, no re-alignment.

BigWig tracks are derived from HSPs by painting each aligned base with its pident value — this is more accurate than the fragment-averaging approach since BLAST computes identity over the true aligned region.

---

## Sensitivity for Short Matches (25-100 bp)

This is the main concern with dropping the fragment approach. Analysis:

| Region size | Identity | Default BLASTN (`word_size=11`) | Tuned (`word_size=7`) |
|-------------|----------|--------------------------------|----------------------|
| 100 bp | >90% | Reliable | Reliable |
| 50 bp | >90% | Usually fine | Reliable |
| 25 bp | >90% | May miss | Usually works |
| 50 bp | ~80% | Can miss | Usually works |
| 25 bp | ~80% | **Likely missed** | Marginal |

**Key insight**: The fragment approach doesn't help much here either. A 50 bp fragment containing a 25 bp match at 80% identity dilutes the signal with 25 bp of unrelated sequence, making it hard to detect regardless.

**Mitigation**: Tuned BLASTN parameters for sensitive/divergent modes:
- `word_size=7` — seeds at every 7-mer exact match
- `-reward 1 -penalty -1` — tolerant scoring for divergent matches
- `-evalue 10` — permissive E-value cutoff
- These are already implemented in `blastn.py` for `--sensitive` / `--divergent` flags

**Practical threshold**: Gene-vs-gene BLAST with tuned parameters reliably catches matches down to ~30 bp at >85% identity. Below that, any method struggles.

---

## What Gets Removed

### Files to Delete

| File | Reason |
|------|--------|
| `crossgene/align.py` | minimap2 wrapper — no longer needed |
| `crossgene/parser.py` | PAF format parser — no longer needed |
| `crossgene/fragment.py` | Rolling-window fragmentation — no longer needed |
| `crossgene/merge.py` | Hit merging — BLAST finds blocks directly |
| `crossgene/scores.py` | Fragment-based per-base scoring — replaced by HSP-based scoring |
| `tests/test_merge.py` | Tests for removed merge module |
| `artifacts/FEATURE_merge_hits.md` | Design doc for removed feature |

### CLI Options to Remove

| Option | Reason |
|--------|--------|
| `--aligner` | Only BLASTN now |
| `--minimap2-preset` | No minimap2 |
| `--fragment-size` | No fragmentation |
| `--step-size` | No fragmentation |
| `--merge` | No merging (BLAST finds blocks directly) |
| `--merge-threshold` | No merging |

### CLI Options to Keep / Modify

| Option | Status | Notes |
|--------|--------|-------|
| `--min-quality` | **Keep** | Filter HSPs by pident |
| `--max-secondary` | **Keep** | Maps to `-max_target_seqs` / `-max_hsps` |
| `--moderate` / `--sensitive` / `--divergent` | **Keep** | Control BLASTN word_size and scoring |
| `--strict` | **Keep** | Adjusts quality thresholds |
| `--min-mapq` | **Keep** | Pseudo-MAPQ derived from bitscore, already implemented in BLAST parser |
| `--blacklist` | **Repurposed** | Hard-mask blacklisted regions with N's in query sequence before BLAST |
| `--bed-all-hits` | **Keep** | Primary vs secondary HSPs |
| `--output-formats` | **Keep** | bigwig, tsv, plot, bed |
| `--bed` / `--bed-color` | **Keep** | Annotation overlay |

### New CLI Options

| Option | Purpose |
|--------|---------|
| `--min-length` | Minimum HSP alignment length to report (filter very short noise hits) |
| `--min-bitscore` | Minimum bitscore to report a hit (replaces MAPQ-based filtering) |
| `--max-evalue` | Maximum E-value to report (BLAST-native quality filter) |

---

## Files to Modify

### `crossgene/blastn.py` — Alignment Engine

**Changes:**
- Move `write_target_fasta()` here (from `align.py`)
- Add `align_genes()` function: takes two `GeneRecord` objects, writes query FASTA, runs BLASTN gene-vs-gene
- Keep existing `BlastParams` and sensitivity presets
- Remove `fragment_size` from `BlastParams` (no longer relevant)
- Remove `align_fragments_blastn()` (fragment-based entry point) — replaced by `align_genes()`

### `crossgene/parser_blast.py` — BLAST Output Parser

**Changes:**
- Add `parse_blast_gene_vs_gene()`: parses BLASTN tabular output where query = full gene (not fragments)
- Coordinate translation: BLAST qstart/qend → genomic query coords, sstart/send → genomic target coords
- No fragment index parsing needed
- Primary/secondary: first HSP is primary, rest are secondary (or use bitscore ranking)
- Remove `parse_blast_tabular()` (fragment-based parser)

### `crossgene/scores.py` — BigWig Scoring (Rewrite)

**Current**: Assumes fixed-size fragments, averages overlapping scores per base.

**New**: Paint per-base scores from HSPs.

```python
def compute_scores_from_hsps(
    hits: list[AlignmentHit],
    query_gene: GeneRecord,
) -> np.ndarray:
    """Compute per-base similarity score from BLAST HSPs.

    For each base in the query gene, assign the identity of the best
    overlapping HSP (or 0 if no HSP covers it). Scale to 0-100.
    """
    gene_len = query_gene.end - query_gene.start
    scores = np.zeros(gene_len, dtype=np.float64)

    for hit in hits:
        start = hit.query_start - query_gene.start
        end = hit.query_end - query_gene.start
        # Best identity wins (max, not average)
        identity_score = hit.identity * 100.0
        mask = scores[start:end] < identity_score
        scores[start:end] = np.where(mask, identity_score, scores[start:end])

    return scores
```

**Why max instead of mean**: With HSPs (unlike overlapping fragments), each base is typically covered by at most one HSP in a given direction. If HSPs overlap, the best identity is the most informative signal. No averaging artifacts.

### `crossgene/cli.py` — CLI Entry Point

**Changes:**
- Remove fragment/merge/minimap2 options and imports
- Simplify `_run_direction()`: extract genes → BLASTN → parse → score → output
- Add new filter options (`--min-length`, `--min-bitscore`, `--max-evalue`)

### `crossgene/models.py` — Data Structures

**Changes:**
- `AlignmentHit`: add `evalue: float` and `bitscore: float` fields (native BLAST metrics)
- `mapq` field: keep for compatibility but set to pseudo-MAPQ derived from bitscore (or remove)

### `crossgene/tsv_writer.py` — TSV Output

**Changes:**
- Add evalue and bitscore columns
- Remove or keep query_coverage (still meaningful: alignment_length / query_length)

---

## New Pipeline Flow

```
Gene A sequence (full)
    │
    ▼
BLASTN: gene A (query) → gene B (subject)
    │
    ▼
Parse HSPs → list[AlignmentHit]
    │
    ├─► Filter by min_quality, min_length, max_evalue
    │
    ├─► compute_scores_from_hsps() → BigWig  (per-base best identity)
    │
    ├─► write_tsv()
    ├─► write_hit_beds()
    └─► circlize plot

(repeat B→A direction)
```

Compare to current:
```
(CURRENT)  fragment → align 100k fragments → parse → score → merge → output
(NEW)      align 1 sequence → parse → filter → score → output
```

---

## Implementation Phases

### Phase 1: Core BLASTN gene-vs-gene pipeline

| Step | Task | Files |
|------|------|-------|
| 1.1 | Add `align_genes()` to `blastn.py` | `blastn.py` |
| 1.2 | Add `parse_blast_gene_vs_gene()` to `parser_blast.py` | `parser_blast.py` |
| 1.3 | Rewrite `scores.py` for HSP-based BigWig scoring | `scores.py` |
| 1.4 | Add `evalue` / `bitscore` fields to `AlignmentHit` | `models.py` |
| 1.5 | Update `tsv_writer.py` for new fields | `tsv_writer.py` |

### Phase 2: CLI simplification

| Step | Task | Files |
|------|------|-------|
| 2.1 | Rewrite `_run_direction()` for gene-vs-gene flow | `cli.py` |
| 2.2 | Remove fragment/merge/minimap2 CLI options | `cli.py` |
| 2.3 | Add `--min-length`, `--min-bitscore`, `--max-evalue` options | `cli.py` |
| 2.4 | Update `--strict` mode defaults for new parameters | `cli.py` |

### Phase 3: Cleanup

| Step | Task | Files |
|------|------|-------|
| 3.1 | Delete `align.py`, `parser.py`, `fragment.py`, `merge.py` | delete files |
| 3.2 | Delete `tests/test_merge.py`, `artifacts/FEATURE_merge_hits.md` | delete files |
| 3.3 | Update/add tests for new pipeline | `tests/` |
| 3.4 | Update `CLAUDE.md` project documentation | `CLAUDE.md` |

---

## Edge Cases

| Case | Handling |
|------|----------|
| **No HSPs found** | Empty hit list → empty TSV, no arcs in plot, flat BigWig. Log warning. |
| **Overlapping HSPs** | BigWig: best identity wins per base. TSV/BED: report all (or filter by primary). |
| **Very long genes (>100 kb)** | Single BLASTN call, no fragmentation overhead. BLAST handles this fine. |
| **Genes on different chromosomes** | Already handled — BLAST doesn't care about chromosomes. |
| **Self-comparison (gene A = gene A)** | BLAST will find a perfect self-hit. Filter or flag it. |
| **Short HSPs (< 25 bp)** | Controlled by `--min-length` filter. Default: 25 bp. |
| **Minus-strand genes** | Sequence already extracted as sense strand. BLAST reports strand of hit. Coordinate translation handles this. |

## Resolved Decisions

1. **`--min-length` default = 25 bp** — matches the smallest sequences of interest.

2. **`--blacklist` repurposed as pre-alignment filter** — hard-mask blacklisted regions in the query gene sequence with N's before writing the BLAST query FASTA. BLASTN cannot seed through N's, so these regions will not produce hits. This is simpler than the current fragment-skipping approach and works at any scale. The `bed_parser.py` module (parse/filter/clip) is reused; only the application changes (mask sequence instead of skip fragments).

3. **`mapq` field kept** — pseudo-MAPQ derived from bitscore is already implemented in `parser_blast.py` and used by downstream filters. No reason to remove it.

## Blacklist Implementation Detail

```python
def mask_sequence(sequence: str, regions: list[BedRegion], gene_start: int) -> str:
    """Replace blacklisted regions with N's in the gene sequence.

    Args:
        sequence: Gene sequence string.
        regions: BedRegions already filtered/clipped to gene bounds.
        gene_start: Genomic start of the gene (for coordinate offset).

    Returns:
        Masked sequence with N's replacing blacklisted bases.
    """
    seq = list(sequence)
    for region in regions:
        local_start = region.start - gene_start
        local_end = region.end - gene_start
        for i in range(max(0, local_start), min(len(seq), local_end)):
            seq[i] = 'N'
    return ''.join(seq)
```

Called in `_run_direction()` before writing the query FASTA, if `--blacklist` is provided.
