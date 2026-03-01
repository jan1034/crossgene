# CrossGene — Architecture

## Design Decisions (from Discussion)

| # | Question | Decision |
|---|----------|----------|
| 1 | Step size | Default = 1 bp, configurable via `--step-size`. Controls visualization density too (no separate subsampling) |
| 2 | Aligner | **minimap2**. Report **all** alignments above threshold, weighted by alignment quality |
| 2b | Multimappers | Keep all hits per fragment; each hit carries its own quality weight |
| 2c | Reverse complement | Yes, detect and flag. Strand used for visualization color |
| 3 | Similarity score | Per-base mappability score, normalized 0–100 |
| 4 | BigWig | Real genomic coordinates from GTF, importable in IGV |
| 5 | TSV columns | Confirmed (see TSV section) |
| 6 | Visualization | Color = strand, intensity = quality. Feature annotations from GTF. Density controlled by step_size |
| 7 | Scope | Gene-level only. No batch mode for now (future consideration) |
| 8 | Direction | **Bidirectional**: A→B and B→A, producing outputs for both directions |
| 9 | Gene name resolution | Search by `gene_name` only. Multiple matches → warn and take first |
| 10 | Strand handling | Use **sense sequence** (revcomp for minus-strand genes). Store and label strand in outputs and visualization |
| 11 | Visualization arcs | **One arc per match** (from A→B direction only, since B→A is largely symmetric). Avoids duplication |
| 12 | minimap2 preset | Auto-select based on `fragment_size`. Primary use case: short fragments (30–100 bp) → `-x sr` |

### Aligner Recommendation: minimap2

minimap2 is preferred over BWA mem for this use case because:
- **Better for divergent sequences**: gene comparisons often involve paralogs or orthologs with significant divergence; minimap2 handles this more robustly
- **Built-in secondary alignments**: `minimap2 -a --secondary=yes -N <max_hits>` natively reports multiple alignment locations
- **Faster indexing**: for single-gene references, indexing is near-instant
- **Flexible presets**: `-x sr` for short fragments (30–100 bp, the primary use case), `-x map-ont` for longer ones
- **No separate index step required**: can align directly against a FASTA reference

**Preset auto-selection** (based on `--fragment-size`):
- ≤200 bp: `-x sr` (short read mode) — **default use case**
- \>200 bp: default preset (no `-x` flag)
- User can override via `--minimap2-preset` if needed

---

## Input Model

### Reference Data (in `references/` directory)
The tool uses a pre-existing genome setup:
- `homo_sapiens.109.mainChr.fa` — whole-genome FASTA (headers: `>chr1`, `>chr2`, ...)
- `homo_sapiens.109.mainChr.fa.fai` — FASTA index
- `homo_sapiens.109.genes.gtf` — gene annotations (Ensembl 109)
- `homo_sapiens.109.chrom.sizes` — chromosome sizes

### User Input
The user provides **gene names** (e.g., `BRCA1`, `BRCA2`). The tool:
1. Looks up the gene in the GTF to get `chr:start-end` and strand
2. Extracts the sequence from the genome FASTA using pysam/samtools
3. Uses the genomic coordinates for all downstream outputs (BigWig, TSV)

This means:
- No separate gene FASTA files needed
- Genomic coordinates are always available (IGV compatibility guaranteed)
- Feature annotations (exons, etc.) are extracted from the same GTF

### Configurable Reference Paths
The paths to genome FASTA, GTF, and chrom.sizes should be configurable (via CLI flags or a config file), defaulting to the `references/` directory.

---

## Bidirectional Comparison

The tool runs the pipeline in **both directions**:

```
Direction 1: Fragment Gene A → Align to Gene B
Direction 2: Fragment Gene B → Align to Gene A
```

This produces:
- **2 BigWig files**: mappability of A relative to B, and B relative to A
- **2 TSV files**: hits from A→B and B→A
- **1 visualization**: single circlize plot showing both directions (or combined arcs)

Regions in A with no match in B show as low-scoring in the A→B BigWig.
Regions in B with no match in A show as low-scoring in the B→A BigWig.
Together this gives a complete picture of shared and unique regions.

---

## Similarity Score — Detailed Explanation

### The Problem
With a rolling window (e.g., 100 bp fragments, step 1 bp), each base in the query gene is covered by up to `fragment_size` overlapping fragments. Each fragment produces 0..N alignments to the target gene. We need to collapse this into a **single per-base score** for the BigWig track.

### Approach: Normalized Mappability Score (0–100)

For each base position `i` in the query gene:

1. Collect all fragments that overlap position `i`
2. For each fragment, take its **best** alignment score (highest identity)
3. Compute the **mean** of these best scores across all overlapping fragments
4. **Normalize to 0–100** scale

```
raw_score(i) = mean( best_identity(f) for f in fragments_covering(i) )
score(i) = raw_score(i) * 100    # identity is 0.0–1.0 → score is 0–100
```

- **100** → perfect match to target gene at this position
- **~80** → good match with some divergence
- **0** → no fragment covering this base aligned at all

Using identity (matches / alignment_length) as the base metric makes the score intuitive and comparable across different fragment sizes.

### What goes into the TSV vs BigWig

| Output | Content |
|--------|---------|
| **BigWig** | Per-base mappability score (0–100). Shows "how well does this region map to the other gene?" |
| **TSV** | Per-fragment, per-alignment detail: all individual hits with coordinates, strand, identity, MAPQ. Enables analysis of multi-mapping and specific region relationships |

---

## Data Flow

```
                         INPUTS
              ┌─────────────────────────┐
              │  Gene A name (e.g. BRCA1) │
              │  Gene B name (e.g. BRCA2) │
              │  Genome FASTA + GTF       │
              └────────────┬──────────────┘
                           │
              ┌────────────▼──────────────┐
              │     Gene Extractor         │
              │  GTF lookup → coordinates  │
              │  FASTA extract → sequence  │
              │  Also extracts exon/feature│
              │  annotations from GTF      │
              └────────────┬──────────────┘
                           │
          ┌────────────────┼────────────────┐
          │ Direction A→B                   │ Direction B→A
          │                                 │
  ┌───────▼──────────┐            ┌────────▼─────────┐
  │ Fragment Gene A   │            │ Fragment Gene B   │
  │ rolling window    │            │ rolling window    │
  └───────┬──────────┘            └────────┬─────────┘
          │                                │
  ┌───────▼──────────┐            ┌────────▼─────────┐
  │ Align → Gene B    │            │ Align → Gene A    │
  │ minimap2          │            │ minimap2           │
  └───────┬──────────┘            └────────┬─────────┘
          │                                │
  ┌───────▼──────────┐            ┌────────▼─────────┐
  │ Parse alignments  │            │ Parse alignments  │
  └───────┬──────────┘            └────────┬─────────┘
          │                                │
  ┌───────▼──────────┐            ┌────────▼─────────┐
  │ BigWig A (0-100)  │            │ BigWig B (0-100)  │
  │ TSV A→B           │            │ TSV B→A           │
  └───────┬──────────┘            └────────┬─────────┘
          │                                │
          └────────────┬───────────────────┘
                       │
              ┌────────▼──────────┐
              │   Visualization    │
              │   pycirclize       │
              │   combined arcs    │
              │   from both dirs   │
              └───────────────────┘
```

---

## Module Boundaries

### 1. `cli.py` — Entry Point
- Parse arguments (click)
- Orchestrate pipeline: extract → fragment → align → parse → output
- Run both directions sequentially

### 2. `gene_extractor.py` — Gene Lookup & Extraction
- Input: gene name, GTF path, genome FASTA path
- Looks up gene in GTF → gets chr, start, end, strand
- Extracts sequence from genome FASTA (using pysam `fetch`)
- Extracts feature annotations (exons, UTRs) from the same GTF for the gene region
- Output: `GeneRecord` with sequence, genomic coordinates, and features
- Searches by `gene_name` attribute only (not gene_id)
- Multiple matches: print warning with all matches, use the first one
- Error handling: gene not found

### 3. `fragment.py` — Fragment Generator
- Input: GeneRecord (sense sequence + coordinates), fragment_size, step_size
- Output: temporary FASTA file with fragments
- Fragment naming convention: `>gene_name:idx` (sequential index, mapped back to genomic coords via offset)
- Uses the **sense sequence** (already reverse-complemented for minus-strand genes by gene_extractor)
- Handles sequences shorter than fragment_size gracefully
- Primary use case: short fragments (30–100 bp)

### 4. `align.py` — Aligner Wrapper
- Input: fragment FASTA, target gene FASTA (extracted), aligner params
- Calls minimap2 via subprocess
- Output: PAF file
- Configurable: max secondary alignments (`-N`), min MAPQ
- **Sensitive mode** (`--sensitive`): uses `-w 5 -k 11` for better detection of short fragment matches (30–50 bp). Default `-x sr` works well for 50–100 bp

**PAF format**: lighter than BAM and sufficient for our needs (coordinates, strand, identity, CIGAR). The BigWig serves the purpose of IGV visualization.

### 5. `parser.py` — Alignment Parser
- Input: PAF file, query gene genomic offset, target gene genomic offset
- Output: list of `AlignmentHit` dataclass objects
- Translates fragment-relative coordinates back to genomic coordinates
- Computes: identity (matches / alignment_length), strand, MAPQ
- Filters by `min_quality` threshold

### 6. `scores.py` — Score Aggregation
- Input: list of `AlignmentHit`, gene length, genomic start offset
- Output: per-base score array (numpy), normalized 0–100
- Aggregation: mean of best-hit identity per overlapping fragment

### 7. `bigwig.py` — BigWig Writer
- Input: per-base score array, chromosome name, genomic start position, chrom.sizes file
- Output: `.bw` file with real genomic coordinates
- Uses pyBigWig
- Chromosome sizes read from the provided `.chrom.sizes` file

### 8. `tsv_writer.py` — TSV Output
- Input: list of `AlignmentHit`, direction label
- Output: TSV file with all hits (one row per alignment)

### 9. `visualize.py` — Circlize Visualization
- Input: AlignmentHits (from A→B direction), gene records (with features + strand info)
- Output: PDF/SVG
- Layout: pycirclize with Gene A right, Gene B left
- **One arc per match** (A→B direction only to avoid near-duplicate arcs)
  - **Color**: alignment strand (e.g., blue = same-sense match, red = antisense/reverse complement match)
  - **Intensity/alpha**: mapping quality (high = opaque, low = transparent)
- Genes displayed in **sense orientation** (5'→3'), labeled with strand: "BRCA1 (+)" or "TP53 (−)"
- Feature track: exons/UTRs from GTF drawn as colored blocks on gene arcs
- Arc density controlled by step_size parameter (step_size=1 → max arcs, step_size=50 → sparser)

### 10. `models.py` — Data Structures

---

## Key Data Structures

```python
@dataclass
class GeneRecord:
    name: str              # gene name (e.g., "BRCA1")
    gene_id: str           # Ensembl ID (e.g., "ENSG00000012048")
    chrom: str             # chromosome (e.g., "chr13")
    start: int             # genomic start (0-based)
    end: int               # genomic end
    strand: str            # '+' or '-' (original strand from GTF)
    sequence: str          # SENSE sequence (revcomp applied if minus strand)
    features: list[GeneFeature]  # exons, UTRs, etc.

@dataclass
class GeneFeature:
    feature_type: str      # "exon", "UTR", "CDS", etc.
    start: int             # genomic start
    end: int               # genomic end
    metadata: dict         # additional GTF attributes

@dataclass
class AlignmentHit:
    # Query gene coordinates (genomic)
    query_chrom: str
    query_start: int       # genomic start of fragment
    query_end: int         # genomic end of fragment

    # Target gene coordinates (genomic)
    target_chrom: str
    target_start: int      # hit start in target (genomic)
    target_end: int        # hit end in target (genomic)

    # Alignment info
    strand: str            # '+' or '-'
    identity: float        # 0.0 - 1.0
    mapq: int              # mapping quality
    cigar: str             # CIGAR string
    alignment_score: int   # AS tag value
    is_primary: bool       # primary vs secondary alignment

    # Metadata
    query_gene: str        # gene name of query
    target_gene: str       # gene name of target
    direction: str         # "A→B" or "B→A"
```

---

## BigWig Coordinate Handling (IGV Compatibility)

Since we extract genes from a genome with known coordinates:

1. Gene lookup in GTF provides `chr`, `start`, `end`
2. Chromosome sizes come from `homo_sapiens.109.chrom.sizes`
3. BigWig entries use absolute genomic positions

**Example**: BRCA1 is on chr17:43044295-43170245.
The BigWig will have entries for `chr17` at positions `43044295..43170245`, each with a mappability score 0–100.

In IGV: load the BigWig → navigate to BRCA1 → see the mappability track aligned with the gene.

**Minus-strand genes**: Fragments are generated from the sense sequence, but coordinates are reverse-mapped to genomic positions. Fragment at sense position 0 → genomic position `end`, sense position N → genomic position `end - N`. BigWig always contains ascending genomic coordinates.

---

## TSV Output Format

```
query_gene  query_chrom  frag_start  frag_end  target_gene  target_chrom  hit_start  hit_end  strand  identity  mapq  alignment_score  is_primary  cigar  direction
BRCA1       chr17        43044295    43044395  BRCA2        chr13         32315480   32315578 +       0.92      40    185              true        100M   A→B
BRCA1       chr17        43044295    43044395  BRCA2        chr13         32350100   32350195 -       0.78      25    140              false       95M5S  A→B
```

- One row per alignment hit (a fragment with 3 hits → 3 rows)
- Coordinates are **genomic**
- Sorted by `frag_start`, then by `alignment_score` descending
- `direction` column distinguishes A→B from B→A hits

---

## Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Gene A | `--gene-a` | required | Gene name (e.g., BRCA1) |
| Gene B | `--gene-b` | required | Gene name (e.g., BRCA2) |
| Fragment size | `--fragment-size` | 50 | Rolling window fragment length (bp). Typical range: 30–100 bp |
| Step size | `--step-size` | 1 | Step between fragments (bp). Also controls visualization density |
| Min quality | `--min-quality` | 30 | Minimum identity (0–100) to report a hit |
| Max secondary | `--max-secondary` | 10 | Max secondary alignments per fragment |
| Genome FASTA | `--genome` | `references/homo_sapiens.109.mainChr.fa` | Path to genome FASTA |
| GTF | `--gtf` | `references/homo_sapiens.109.genes.gtf` | Path to gene annotation GTF |
| Chrom sizes | `--chrom-sizes` | `references/homo_sapiens.109.chrom.sizes` | Chromosome sizes file |
| Output dir | `--outdir` | `.` | Output directory |
| Formats | `--output-formats` | all | Comma-separated: bigwig,tsv,plot |
| minimap2 preset | `--minimap2-preset` | auto | Override minimap2 preset (auto-selects `-x sr` for ≤200 bp) |
| Sensitive mode | `--sensitive` | false | Tune minimap2 for short fragments (≤50 bp): lower minimizer window (`-w 5 -k 11`) for better sensitivity at cost of speed |

---

## Output Files

For a comparison of BRCA1 vs BRCA2:
```
outdir/
├── BRCA1_vs_BRCA2.mappability.bw      # BigWig: BRCA1 mappability to BRCA2
├── BRCA2_vs_BRCA1.mappability.bw      # BigWig: BRCA2 mappability to BRCA1
├── BRCA1_vs_BRCA2.hits.tsv            # TSV: all A→B hits
├── BRCA2_vs_BRCA1.hits.tsv            # TSV: all B→A hits
└── BRCA1_vs_BRCA2.circlize.pdf        # Visualization: combined plot
```

---

## File Structure

```
crossgene/
├── crossgene/
│   ├── __init__.py
│   ├── cli.py              # Click CLI entry point
│   ├── gene_extractor.py   # GTF lookup + FASTA extraction
│   ├── fragment.py          # Rolling window fragmentation
│   ├── align.py             # minimap2 wrapper
│   ├── parser.py            # PAF parsing
│   ├── scores.py            # Per-base score aggregation
│   ├── bigwig.py            # BigWig writer (pyBigWig)
│   ├── tsv_writer.py        # TSV output
│   ├── visualize.py         # Pycirclize visualization
│   └── models.py            # GeneRecord, AlignmentHit dataclasses
├── references/
│   ├── homo_sapiens.109.mainChr.fa
│   ├── homo_sapiens.109.mainChr.fa.fai
│   ├── homo_sapiens.109.genes.gtf
│   └── homo_sapiens.109.chrom.sizes
├── tests/
│   ├── test_fragment.py
│   ├── test_parser.py
│   ├── test_scores.py
│   └── test_data/
├── ANALYSIS.md
├── ARCHITECTURE.md
├── IMPLEMENTATION.md
└── pyproject.toml
```

---

## Technology Stack

| Component | Tool/Library |
|-----------|-------------|
| Language | Python 3.10+ |
| CLI | click |
| Gene extraction | pysam (FASTA fetch), GTF parsing (gtfparse or custom) |
| Fragmentation | Plain Python |
| Alignment | minimap2 (subprocess) |
| PAF parsing | Custom (tab-delimited, straightforward) |
| Score arrays | numpy |
| BigWig writing | pyBigWig |
| TSV | csv (stdlib) |
| Visualization | pycirclize, matplotlib |

---

## Resolved Questions

| # | Question | Resolution |
|---|----------|------------|
| Q1 | Gene name resolution | Search by `gene_name` only. Multiple matches → warn, take first |
| Q2 | Gene strand handling | Use sense sequence (revcomp for minus strand). Store original strand, label in plot |
| Q3 | Visualization arcs | One arc per match (A→B only). Avoids duplication since B→A is largely symmetric |
| Q4 | minimap2 preset | Auto-select based on fragment_size (≤200 bp → `-x sr`). User can override via `--minimap2-preset` |
| Q5 | Minus-strand coordinates | BigWig/TSV always use **genomic positions** (ascending). Internal reverse mapping for minus-strand genes |
| Q6 | Visualization orientation | **Sense orientation** (5'→3') with strand label (+/−) |
| Q7 | minimap2 sensitivity | `--sensitive` flag for short fragments (lowers minimizer window size) |

---

## All questions resolved — no open questions remain.
