# Compare Genes — Implementation Plan

## Design Decisions (from Implementation Discussion)

| # | Question | Decision |
|---|----------|----------|
| I1 | GTF parsing | Use **gtfparse** library (returns pandas DataFrame) |
| I2 | Temporary files | **tempfile** module with auto-cleanup after pipeline completes |
| I3 | Testing strategy | Small synthetic sequences for unit tests, mock subprocess for align module |
| I4 | Plot format | **PDF** only (no format flag for now) |
| I5 | Arc limit | **Auto-subsample**: if arc count > threshold (5000), subsample keeping top-scoring arcs |
| I6 | Logging | **Python logging** module with `--verbose` flag (default: INFO, verbose: DEBUG) |
| I7 | Packaging | **pip-installable** CLI via console_scripts entry point (`compare-genes`) |

---

## Steps

### Step 1 — Project Scaffolding

Create the package skeleton and build configuration.

**Files to create:**
- `pyproject.toml` — package metadata, dependencies, console_scripts entry point
- `compare_genes/__init__.py` — empty (or version only)
- `compare_genes/models.py` — dataclasses (GeneRecord, GeneFeature, AlignmentHit)

**pyproject.toml details:**
- Build backend: `hatchling` (or `setuptools`)
- Python: `>=3.10`
- Dependencies: `click`, `pysam`, `numpy`, `pyBigWig`, `pycirclize`, `gtfparse`, `matplotlib`
- Entry point: `compare-genes = compare_genes.cli:main`

**Checkpoint:** `pip install -e .` succeeds, `compare-genes --help` shows placeholder help text.

---

### Step 2 — Data Models (`models.py`)

Implement the three core dataclasses from ARCHITECTURE.md:
- `GeneRecord` (name, gene_id, chrom, start, end, strand, sequence, features)
- `GeneFeature` (feature_type, start, end, metadata)
- `AlignmentHit` (query coords, target coords, alignment info, metadata)

No external dependencies for this module.

**Checkpoint:** Import and instantiate all dataclasses in a Python shell.

---

### Step 3 — Gene Extractor (`gene_extractor.py`)

Implement gene lookup and sequence extraction.

**Functions:**
- `lookup_gene(gene_name: str, gtf_path: str) -> GeneRecord` (without sequence)
  - Parse GTF with gtfparse, filter by `gene_name` attribute
  - Extract gene-level row for coordinates (chr, start, end, strand)
  - Extract feature rows (exon, CDS, UTR) as `GeneFeature` list
  - Multiple matches: log warning, take first
  - Gene not found: raise descriptive error
- `extract_sequence(gene: GeneRecord, genome_path: str) -> GeneRecord`
  - Open FASTA with `pysam.FastaFile`
  - Fetch sequence for `chrom:start-end`
  - If minus strand: reverse-complement the sequence
  - Return updated GeneRecord with sequence populated

**Tests (`tests/test_gene_extractor.py`):**
- Create a tiny GTF and FASTA in test fixtures
- Test lookup: correct gene found, multiple matches warning, gene not found error
- Test extraction: correct sequence, reverse-complement for minus strand

**Checkpoint:** Extract a real gene (e.g., a small gene) from the reference files and print its length + first 50bp.

---

### Step 4 — Fragment Generator (`fragment.py`)

Implement rolling-window fragmentation.

**Functions:**
- `generate_fragments(gene: GeneRecord, fragment_size: int, step_size: int) -> Path`
  - Create temp FASTA file (tempfile.NamedTemporaryFile, delete=False)
  - Write fragments as `>gene_name:idx` with idx = 0, 1, 2, ...
  - Handle edge case: sequence shorter than fragment_size (single fragment = full sequence)
  - Handle edge case: last fragment may be shorter (include it if >= fragment_size / 2, skip otherwise)
  - Return path to temp FASTA
- `fragment_index_to_genomic(gene: GeneRecord, idx: int, fragment_size: int, step_size: int) -> tuple[int, int]`
  - Convert fragment index to genomic coordinates
  - Account for strand direction (minus-strand: map sense position to descending genomic coords)

**Tests (`tests/test_fragment.py`):**
- Known sequence → verify fragment count, content, naming
- Edge cases: short sequence, step > 1, minus-strand coordinate mapping

**Checkpoint:** Fragment a 500bp test sequence with fragment_size=50, step=10, verify output FASTA has correct number of entries.

---

### Step 5 — Aligner Wrapper (`align.py`)

Implement minimap2 subprocess call.

**Functions:**
- `align_fragments(fragment_fasta: Path, target_fasta: Path, params: AlignParams) -> Path`
  - Write target gene sequence to a temp FASTA (for minimap2 reference)
  - Build minimap2 command:
    - Base: `minimap2 -c --eqx --secondary=yes -N {max_secondary}`
    - Preset: auto-select based on fragment_size, or user override
    - Sensitive mode: add `-w 5 -k 11`
    - Output: PAF format (default, no `-a` flag)
  - Run via `subprocess.run`, capture stderr for logging
  - Check return code, raise on failure
  - Return path to PAF output file
- `AlignParams` dataclass: fragment_size, max_secondary, minimap2_preset, sensitive

**Validation:**
- Check minimap2 is installed (`shutil.which('minimap2')`) at startup, clear error if missing

**Tests (`tests/test_align.py`):**
- Mock subprocess.run, verify command construction
- Test preset auto-selection logic
- Test sensitive mode flags

**Checkpoint:** Align synthetic fragments against a synthetic target, verify PAF output is non-empty.

---

### Step 6 — PAF Parser (`parser.py`)

Parse minimap2 PAF output and translate coordinates.

**Functions:**
- `parse_paf(paf_path: Path, query_gene: GeneRecord, target_gene: GeneRecord, min_quality: int, direction: str) -> list[AlignmentHit]`
  - Read PAF line by line (tab-delimited, 12+ columns)
  - Extract: query name, query length, query start/end, strand, target name, target length, target start/end, matches, alignment_length, mapq
  - Compute identity = matches / alignment_length
  - Filter by min_quality (identity * 100 >= min_quality)
  - Translate query fragment-local coords to genomic coords using fragment index from name
  - Translate target-local coords to genomic coords using target gene offset
  - Extract CIGAR from optional fields (cg:Z: tag)
  - Extract alignment score from optional fields (AS:i: tag)
  - Determine is_primary from tp:A: tag

**PAF columns (0-indexed):**
0: query_name, 1: query_length, 2: query_start, 3: query_end, 4: strand,
5: target_name, 6: target_length, 7: target_start, 8: target_end,
9: num_matches, 10: alignment_block_length, 11: mapq

**Tests (`tests/test_parser.py`):**
- Hand-crafted PAF lines → verify parsed AlignmentHit fields
- Test coordinate translation (plus and minus strand genes)
- Test filtering by quality threshold

**Checkpoint:** Parse a real PAF from step 5, print hit count and first few hits.

---

### Step 7 — Score Aggregation (`scores.py`)

Compute per-base mappability scores.

**Functions:**
- `compute_scores(hits: list[AlignmentHit], gene: GeneRecord, fragment_size: int, step_size: int) -> np.ndarray`
  - Create array of length `gene.end - gene.start`, initialized to 0
  - Group hits by fragment index (derived from query_start)
  - For each fragment: take best identity (highest)
  - For each base covered by that fragment: accumulate best identity and count
  - Final score per base = (sum of best identities / count) * 100
  - Bases with no covering fragments get score 0

**Tests (`tests/test_scores.py`):**
- Synthetic hits with known positions → verify per-base scores
- Test overlapping fragments average correctly
- Test uncovered bases = 0

**Checkpoint:** Compute scores for a small test case, verify array values manually.

---

### Step 8 — BigWig Writer (`bigwig.py`)

Write per-base scores as a BigWig file.

**Functions:**
- `write_bigwig(scores: np.ndarray, chrom: str, start: int, chrom_sizes_path: str, output_path: str) -> None`
  - Read chrom.sizes file into dict
  - Open BigWig for writing with pyBigWig
  - Add chromosome headers
  - Write scores as intervals: each base position `start + i` gets `scores[i]`
  - Use `bw.addEntries` with arrays for efficiency (not one entry at a time)

**Tests (`tests/test_bigwig.py`):**
- Write a small BigWig, read it back with pyBigWig, verify values

**Checkpoint:** Write a BigWig for a test gene, open in IGV or verify with pyBigWig.

---

### Step 9 — TSV Writer (`tsv_writer.py`)

Write alignment hits as a TSV file.

**Functions:**
- `write_tsv(hits: list[AlignmentHit], output_path: str) -> None`
  - Write header row (15 columns from ARCHITECTURE.md)
  - Sort hits by frag_start, then alignment_score descending
  - Write one row per AlignmentHit
  - Use csv.writer with `\t` delimiter

**Tests (`tests/test_tsv_writer.py`):**
- Write TSV from synthetic hits, read back and verify columns/sorting

**Checkpoint:** Write a TSV, inspect output manually.

---

### Step 10 — Visualization (`visualize.py`)

Create circlize plot with pycirclize.

**Functions:**
- `create_circlize_plot(hits_ab: list[AlignmentHit], gene_a: GeneRecord, gene_b: GeneRecord, output_path: str, max_arcs: int = 5000) -> None`
  - Set up pycirclize Circos with two sectors (Gene A right, Gene B left)
  - Draw gene tracks with feature annotations (exons as colored blocks)
  - Label genes with strand info: "BRCA1 (+)" or "TP53 (−)"
  - Draw arcs from A→B hits:
    - Color by strand (blue = same-sense, red = antisense)
    - Alpha by mapq (normalized to 0.2–1.0 range)
  - **Auto-subsample**: if len(hits_ab) > max_arcs, keep top `max_arcs` by alignment_score
  - Save as PDF

**Tests (`tests/test_visualize.py`):**
- Smoke test: create plot from synthetic data, verify PDF file is created and non-empty

**Checkpoint:** Generate a plot from real data, visually inspect.

---

### Step 11 — CLI Entry Point (`cli.py`)

Wire everything together.

**Functions:**
- `main()` — Click command with all parameters from ARCHITECTURE.md
- Pipeline orchestration:
  1. Validate inputs (minimap2 installed, files exist)
  2. Set up logging (INFO default, DEBUG with `--verbose`)
  3. Extract Gene A and Gene B (gene_extractor)
  4. **Direction A→B:**
     a. Fragment Gene A
     b. Write Gene B sequence to temp FASTA (for minimap2 reference)
     c. Align fragments → Gene B
     d. Parse PAF → AlignmentHit list
     e. Compute scores
     f. Write BigWig (if format selected)
     g. Write TSV (if format selected)
  5. **Direction B→A:** repeat steps 4a–4g swapped
  6. **Visualization** (if format selected): create circlize plot from A→B hits
  7. Clean up temp files
  8. Log summary: output file paths, hit counts, runtime

**Click parameters:**
- All from the ARCHITECTURE.md parameters table
- Add: `--verbose` / `-v` flag for debug logging

**Error handling:**
- Gene not found → clear message with suggestions (check gene name spelling)
- minimap2 not installed → clear message with install instructions
- No alignments found → warn but still produce empty outputs
- Output dir doesn't exist → create it

**Tests (`tests/test_cli.py`):**
- Test CLI argument parsing (click.testing.CliRunner)
- Integration smoke test with tiny synthetic data

**Checkpoint:** Full end-to-end run with two real genes. Verify all 5 output files are produced.

---

### Step 12 — End-to-End Testing & Polish

- Run with a real gene pair (e.g., BRCA1 vs BRCA2) and inspect all outputs
- Verify BigWig loads correctly in IGV
- Verify TSV is well-formed
- Verify circlize plot is visually correct
- Profile performance for large genes, optimize if needed (e.g., batch pyBigWig writes)
- Check edge cases: very small genes, genes on different strands, same gene vs itself
