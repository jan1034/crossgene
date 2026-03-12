# Plan: Codebase Cleanup (post revcomp removal)

Apply after `PLAN_remove_revcomp.md` is complete. Grouped by file.

**Baseline:** 95/96 tests pass. 1 known failure (`test_max_secondary`) fixed in step 3.

---

## 1. `models.py` — Remove dead fields

### 1a. Delete `cigar` from `AlignmentHit`

`models.py:65`: `cigar: str` is always `""` (set in `parser_blast.py:182`), never read anywhere.

```python
# DELETE from AlignmentHit:
cigar: str  # CIGAR string
```

Remove the corresponding `cigar=""` in `parser_blast.py:182`.

**Tests to update** (all construct `AlignmentHit` with `cigar=`):
- `test_bed_writer.py:29` — remove `cigar="50M"` from `_make_hit`
- `test_tsv_writer.py:13` — remove `cigar=""` from `_make_hit` defaults
- `test_scores.py:21` — remove `cigar=""` from `_make_hit` defaults
- `test_visualize.py:31` — remove `cigar="50M"` from `_make_hit`

### 1b. Delete `alignment_score`, keep `bitscore`

`models.py:66`: `alignment_score: int` is just `int(bitscore)` (set in `parser_blast.py:183`). The float `bitscore` field already exists and is more precise.

```python
# DELETE from AlignmentHit:
alignment_score: int  # AS tag value or bitscore
```

**Source references to update:**
- `tsv_writer.py:13` — remove `"alignment_score"` from `COLUMNS`
- `tsv_writer.py:44` — remove `h.alignment_score` from the row
- `visualize.py:79` — `_subsample_hits` sorts by `h.alignment_score` → change to `h.bitscore`
- `visualize.py:72` — update docstring from "alignment_score" to "bitscore"
- `parser_blast.py:183` — remove `alignment_score=int(bitscore)`

**Tests to update** (all construct `AlignmentHit` with `alignment_score=`):
- `test_bed_writer.py:29` — remove `alignment_score=100` from `_make_hit`
- `test_tsv_writer.py:9` — `_make_hit` takes `alignment_score` as positional arg → remove it, update signature
- `test_tsv_writer.py:35` — `len(header) == 17` → change to `16`
- `test_tsv_writer.py:44-57` — column index assertions shift left by 1 after removed column (row[14] was `alignment_score`)
- `test_scores.py:21` — remove `alignment_score=90` from `_make_hit` defaults
- `test_visualize.py:27` — `_make_hit` has `score` param → `alignment_score=score` → remove both

### 1c. Delete `metadata` from `GeneFeature`

`models.py:27`: `metadata: dict` is always `{}` (set in `gene_extractor.py:69,155`), never read.

```python
# DELETE from GeneFeature:
metadata: dict = field(default_factory=dict)
```

**Source references to update:**
- `gene_extractor.py:69` — remove `metadata={}` from `GeneFeature(...)` in `lookup_gene`
- `gene_extractor.py:155` — remove `metadata={}` from `GeneFeature(...)` in `load_features`

**Tests to update:**
- `test_visualize.py:74` — `GeneFeature("exon", 100, 200, {})` → remove positional `{}`
- `test_visualize.py:99` — `GeneFeature("exon", 100, 200, {})` → remove positional `{}`

---

## 2. `blastn.py` — Use `dataclasses.replace` for masked gene

`blastn.py:192-199`: Manual field-by-field copy is fragile.

```python
# BEFORE:
masked_gene = GeneRecord(
    name=query_gene.name, gene_id=query_gene.gene_id,
    chrom=query_gene.chrom, start=query_gene.start, end=query_gene.end,
    strand=query_gene.strand, sequence=masked_seq,
    features=query_gene.features,
    gene_body_start=query_gene.gene_body_start,
    gene_body_end=query_gene.gene_body_end,
)

# AFTER:
from dataclasses import replace
masked_gene = replace(query_gene, sequence=masked_seq)
```

---

## 3. `blastn.py` — Fix `max_target_seqs` hardcoded to 1

`blastn.py:87`: Currently hardcoded. Should use the param.

```python
# BEFORE:
"-max_target_seqs", "1",

# AFTER:
"-max_target_seqs", str(params.max_secondary),
```

This fixes `test_blastn.py::TestBuildBlastnCommand::test_max_secondary`.

---

## 4. `bed_writer.py` — Deduplicate query/target BED blocks

`bed_writer.py:51-82`: Two near-identical blocks. Extract a helper.

Use a simple helper that takes concrete values rather than attribute names:

```python
def _write_bed_file(
    path: str, track_name: str, track_desc: str,
    sorted_hits: list[AlignmentHit],
    get_coords: Callable[[AlignmentHit], tuple[str, int, int]],
    primary_prefix: str, secondary_prefix: str,
    primary_num: dict[int, int], secondary_num: dict[int, int],
) -> None:
    with open(path, "w") as f:
        f.write(f'track name="{track_name}" description="{track_desc}" itemRgb=on\n')
        for h in sorted_hits:
            pn = primary_num[id(h)]
            sn = secondary_num[id(h)]
            name = f"{primary_prefix}{pn}-{secondary_prefix}{sn}"
            score = min(int(h.identity * 1000), 1000)
            rgb = COLOR_SENSE if h.strand == "+" else COLOR_ANTISENSE
            chrom, start, end = get_coords(h)
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{h.strand}\t{start}\t{end}\t{rgb}\n")
```

Then call it twice in `write_hit_beds`:

```python
_write_bed_file(
    q_bed_path, ..., query_sorted,
    get_coords=lambda h: (h.query_chrom, h.query_start, h.query_end),
    ...
)
_write_bed_file(
    t_bed_path, ..., target_sorted,
    get_coords=lambda h: (h.target_chrom, h.target_start, h.target_end),
    ...
)
```

Note: the `id()` usage is fragile but functional here because all references are within the same function scope (no GC possible). Not worth changing unless the structure is already being refactored.

---

## 5. ~~`bed_writer.py` — Drop unused return value~~ — SKIP

`bed_writer.py:87` returns `(q_bed_path, t_bed_path)`. `cli.py` discards it but
all tests in `test_bed_writer.py` capture and use the return value (`q_bed, t_bed = write_hit_beds(...)`).
Removing the return creates test churn for zero benefit. **Keep the return as-is.**

---

## 6. `parser_blast.py` — Single-pass BLAST parsing

`parser_blast.py:84-96`: Two passes over the lines list — first for max bitscore, then to build hits.

Merge into one pass: build hits first, then normalize MAPQ.

```python
# AFTER:
hits: list[AlignmentHit] = []
max_bs = 0.0

for line in lines:
    # ... parse fields, apply filters (except min_mapq) ...
    if bitscore > max_bs:
        max_bs = bitscore
    hits.append(AlignmentHit(..., mapq=0))  # placeholder

# Normalize MAPQ
if max_bs == 0:
    max_bs = 1.0
for h in hits:
    h.mapq = min(60, int(h.bitscore / max_bs * 60))

# Apply min_mapq filter
hits = [h for h in hits if h.mapq >= min_mapq]
```

This eliminates the line-storage list and second iteration. Note: `AlignmentHit` is a dataclass (not frozen), so field assignment works.

---

## 7. `cli.py` — Bundle filter thresholds into a dataclass

`cli.py:45-50`: `_run_direction` has 14 parameters, most are filter thresholds passed through.

```python
# In models.py:
@dataclass
class FilterParams:
    min_quality: float = 30
    min_mapq: int = 0
    min_length: int = 25
    min_bitscore: float = 0.0
    max_evalue: float = float("inf")
```

Then `_run_direction` becomes:

```python
def _run_direction(
    query_gene, target_gene, blast_params, filters: FilterParams,
    chrom_sizes, outdir, formats, direction_label,
    blacklist_regions=None, bed_primary_only=False,
) -> tuple[list, list[Path]]:
```

And in `parse_blast_gene_vs_gene`:

```python
hits = parse_blast_gene_vs_gene(
    tsv_path, query_gene, target_gene,
    filters.min_quality, direction_label,
    min_mapq=filters.min_mapq,
    min_length=filters.min_length,
    min_bitscore=filters.min_bitscore,
    max_evalue=filters.max_evalue,
)
```

Update call sites in `cli.py:256-272` to construct `FilterParams` once and pass it.

---

## 8. `cli.py` — Move `bed_writer` import to top level

`cli.py:97`: Late import inside `_run_direction`.

```python
# BEFORE (inside function):
if "bed" in formats:
    from crossgene.bed_writer import write_hit_beds

# AFTER (top of file, with other imports):
from crossgene.bed_writer import write_hit_beds
```

---

## 9. `bigwig.py` — Accept chrom_sizes dict instead of path

`bigwig.py:48`: `write_bigwig` reads the chrom.sizes file on every call (called twice — once per direction).

```python
# BEFORE:
def write_bigwig(scores, chrom, start, chrom_sizes_path: str, output_path):
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)

# AFTER:
def write_bigwig(scores, chrom, start, chrom_sizes: dict[str, int], output_path):
    ...
```

In `cli.py`, parse once before the direction calls:

```python
from crossgene.bigwig import read_chrom_sizes, write_bigwig
chrom_sizes_dict = read_chrom_sizes(chrom_sizes)
```

Then pass `chrom_sizes_dict` to `_run_direction` and through to `write_bigwig`.

**Tests to update** (`test_bigwig.py` passes file paths to `write_bigwig`):
- `test_bigwig.py:42` — `test_write_and_read_back`: pass dict instead of path
- `test_bigwig.py:53` — `test_larger_array`: pass dict instead of path
- `test_bigwig.py:66` — `test_unknown_chrom_raises`: pass dict instead of path
- `test_bigwig.py:72` — `test_empty_scores`: pass dict instead of path
- Fixtures `chrom_sizes_file` / `chrom_sizes_with_header` are still useful for testing `read_chrom_sizes` itself. Add a `chrom_sizes_dict` fixture for `write_bigwig` tests.

---

## 10. `visualize.py` — Replace defensive bounds clamping with warning + skip

`visualize.py:200-204`: Clamps arc coords to sector bounds. BLAST output is bounded by gene length, so this should never trigger. If it does, it masks an upstream bug.

```python
# BEFORE:
a_local_start = max(0, min(a_local_start, gene_a_len))
a_local_end = max(0, min(a_local_end, gene_a_len))
b_local_start = max(0, min(b_local_start, gene_b_len))
b_local_end = max(0, min(b_local_end, gene_b_len))

# AFTER:
if not (0 <= a_local_start <= a_local_end <= gene_a_len):
    logger.warning("Hit coords out of bounds for %s: [%d, %d) in gene len %d — skipping",
                   gene_a.name, a_local_start, a_local_end, gene_a_len)
    continue
if not (0 <= b_local_start <= b_local_end <= gene_b_len):
    logger.warning("Hit coords out of bounds for %s: [%d, %d) in gene len %d — skipping",
                   gene_b.name, b_local_start, b_local_end, gene_b_len)
    continue
```

This makes out-of-bounds visible instead of silently hiding it.

---

## 11. `gene_extractor.py` — Use `dataclasses.replace` in `load_features`

`gene_extractor.py:165-176`: Manual field-by-field copy of `GeneRecord`, same pattern as step 2.

```python
# BEFORE:
return GeneRecord(
    name=gene.name, gene_id=gene.gene_id, chrom=gene.chrom,
    start=gene.start, end=gene.end, strand=gene.strand,
    sequence=gene.sequence, features=features,
    gene_body_start=gene.gene_body_start, gene_body_end=gene.gene_body_end,
)

# AFTER:
from dataclasses import replace
return replace(gene, features=features)
```

---

## 12. `gene_extractor.py` — Use `dataclasses.replace` in `extract_sequence`

`gene_extractor.py:214-225`: Manual field-by-field copy of `GeneRecord`, same pattern.

```python
# BEFORE:
return GeneRecord(
    name=gene.name, gene_id=gene.gene_id, chrom=gene.chrom,
    start=flanked_start, end=flanked_end, strand=gene.strand,
    sequence=seq, features=gene.features,
    gene_body_start=gene_body_start, gene_body_end=gene_body_end,
)

# AFTER:
from dataclasses import replace
return replace(gene, start=flanked_start, end=flanked_end,
               sequence=seq, gene_body_start=gene_body_start,
               gene_body_end=gene_body_end)
```

---

## 13. `gene_extractor.py` — Remove dead feature extraction in `lookup_gene`

`gene_extractor.py:56-71`: `lookup_gene` extracts sub-gene features (exon, CDS, UTR) from the genes GTF. But when visualization is requested, `load_features` in `cli.py:216-217` re-extracts features from the annotation GTF and completely replaces the features list. When visualization is NOT requested, features are never used.

This code is dead weight — it parses rows and builds `GeneFeature` objects that are always either overwritten or ignored.

```python
# DELETE from lookup_gene (lines 56-71):
    # Extract sub-gene features (exon, CDS, UTR, etc.) if available
    features = []
    gene_id = str(gene_row["gene_id"])
    sub_features = df[
        (df["gene_id"] == gene_id) & (df["feature"] != "gene")
    ]
    for _, feat_row in sub_features.iterrows():
        feature_type = str(feat_row["feature"])
        if feature_type in ("exon", "CDS", "five_prime_utr", "three_prime_utr"):
            features.append(
                GeneFeature(
                    feature_type=feature_type,
                    start=int(feat_row["start"]) - 1,
                    end=int(feat_row["end"]),
                    metadata={},
                )
            )

# And simplify the return:
    gene_id = str(gene_row["gene_id"])
    return GeneRecord(
        name=gene_name,
        gene_id=gene_id,
        chrom=str(gene_row["seqname"]),
        start=start_0based,
        end=end_0based,
        strand=str(gene_row["strand"]),
        sequence="",
    )
```

This also eliminates the `GeneFeature` import from `lookup_gene`'s scope if it's only needed by `load_features`.

Note: the genes GTF (`homo_sapiens.109.genes.gtf`) may not even contain sub-gene features — it's the *annotation* GTF (`mainChr.gtf`) that does. This further confirms the code is dead.

---

## 14. `cli.py` — Minor: `sum()` over tuple instead of list

`cli.py:227`: `sum([moderate, sensitive, divergent])` allocates an unnecessary list.

```python
# BEFORE:
mode_count = sum([moderate, sensitive, divergent])

# AFTER:
mode_count = sum((moderate, sensitive, divergent))
```

Trivial, can be done in any commit that touches `cli.py`.

---

## Implementation order

| Step | Change | Risk | Files touched |
|------|--------|------|---------------|
| 3 | Fix `max_target_seqs` bug | None — fixes existing failing test | `blastn.py` |
| 2 | `dataclasses.replace` for masked gene | None — same behavior | `blastn.py` |
| 11 | `dataclasses.replace` in `load_features` | None — same behavior | `gene_extractor.py` |
| 12 | `dataclasses.replace` in `extract_sequence` | None — same behavior | `gene_extractor.py` |
| 8 | Move `bed_writer` import to top | None — style only | `cli.py` |
| 14 | `sum()` tuple fix | None — trivial | `cli.py` |
| 1a | Remove `cigar` | None — dead code | `models.py`, `parser_blast.py`, 4 test files |
| 1c | Remove `metadata` | None — dead code | `models.py`, `gene_extractor.py`, `test_visualize.py` |
| 13 | Remove dead feature extraction in `lookup_gene` | Low — verify genes GTF has no sub-gene features | `gene_extractor.py` |
| 1b | Remove `alignment_score` | Low — update TSV columns, visualize sort, 4 test files | `models.py`, `parser_blast.py`, `tsv_writer.py`, `visualize.py`, 4 test files |
| 4 | Deduplicate BED writing | Low — extract helper | `bed_writer.py` |
| 9 | Pass chrom_sizes as dict | Low — change function signatures | `bigwig.py`, `cli.py`, `test_bigwig.py` |
| 6 | Single-pass BLAST parsing | Medium — restructure parse loop | `parser_blast.py` |
| 7 | `FilterParams` dataclass | Medium — touches cli.py and parser_blast.py | `models.py`, `cli.py` |
| 10 | Bounds clamping → warning + skip | Low — changes error handling | `visualize.py` |

Do the no-risk items first, then work down. Steps 2+11+12 can be one commit (all `dataclasses.replace`). Steps 1a+1c can be one commit (dead field removal). Each remaining step should be a separate commit to keep changes reviewable.

---

## Suggested commit grouping

1. **Fix `max_target_seqs` bug** (step 3)
2. **Use `dataclasses.replace` everywhere** (steps 2, 11, 12)
3. **Remove dead fields: `cigar`, `metadata`** (steps 1a, 1c)
4. **Remove dead feature extraction in `lookup_gene`** (step 13)
5. **Remove `alignment_score` field** (step 1b) — larger change, separate commit
6. **Import and style cleanups** (steps 8, 14)
7. **Deduplicate BED writing** (step 4)
8. **Pass chrom_sizes as dict** (step 9)
9. **Single-pass BLAST parsing** (step 6)
10. **Bundle filter thresholds into `FilterParams`** (step 7)
11. **Bounds clamping → warning + skip** (step 10)

---

## Files NOT changed

- `scores.py` — clean, no issues
- `bed_parser.py` — clean, no issues
- `bed_writer.py` return value — kept as-is (tests rely on it, removing adds churn for no benefit)
