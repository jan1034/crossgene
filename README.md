# CrossGene

A Python CLI tool for comparing two gene sequences using direct BLASTN alignment (gene-vs-gene). It aligns one gene against another using [BLASTN](https://blast.ncbi.nlm.nih.gov/), identifies homologous regions (HSPs), and produces similarity profiles, detailed hit tables, BED files, and circular visualizations. The comparison runs bidirectionally (A→B and B→A).

## Outputs

For a comparison of e.g. BRCA1 vs BRCA2, the tool produces:

| File | Format | Description |
|------|--------|-------------|
| `BRCA1_vs_BRCA2.mappability.bw` | BigWig | Per-base similarity score (0–100) of BRCA1 relative to BRCA2 |
| `BRCA2_vs_BRCA1.mappability.bw` | BigWig | Per-base similarity score (0–100) of BRCA2 relative to BRCA1 |
| `BRCA1_vs_BRCA2.hits.tsv` | TSV | All individual HSPs from A→B |
| `BRCA2_vs_BRCA1.hits.tsv` | TSV | All individual HSPs from B→A |
| `BRCA1_hits_AtoB.bed` | BED9 | Query-side alignment hits (A→B) |
| `BRCA2_hits_AtoB.bed` | BED9 | Target-side alignment hits (A→B) |
| `BRCA2_hits_BtoA.bed` | BED9 | Query-side alignment hits (B→A) |
| `BRCA1_hits_BtoA.bed` | BED9 | Target-side alignment hits (B→A) |
| `BRCA1_vs_BRCA2.circlize.pdf` | PDF | Circular plot with arcs connecting matching regions |

BigWig files use real genomic coordinates and can be loaded directly into IGV. BED9 files include color coding (steelblue for same-sense, firebrick for antisense) and cross-reference numbering between query and target hits.

## Installation

### Prerequisites

- Python 3.10+
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) on PATH (`blastn` and `makeblastdb`)
- Reference data (genome FASTA, GTF annotations, chromosome sizes)

### Install

```bash
conda create -n crossgene python=3.11
conda activate crossgene
conda install -c bioconda blast

git clone <repository-url>
cd crossgene
pip install -e .
```

Verify the installation:

```bash
crossgene --help
```

## Reference Data

The tool requires a genome FASTA, gene annotation GTF, and chromosome sizes file. By default it looks for these in a `references/` directory:

```
references/
├── homo_sapiens.109.mainChr.fa       # Genome FASTA
├── homo_sapiens.109.mainChr.fa.fai   # FASTA index
├── homo_sapiens.109.genes.gtf        # Gene annotations (Ensembl 109, for gene lookup)
├── homo_sapiens.109.mainChr.gtf      # Annotation GTF (for sub-gene features: exons, CDS, UTRs)
└── homo_sapiens.109.chrom.sizes      # Chromosome sizes
```

These paths are configurable via CLI flags (see below).

## Usage

### Basic usage

```bash
crossgene --gene-a BRCA1 --gene-b BRCA2
```

### Divergent paralogs (high sensitivity)

Use low sensitivity and stringency levels for divergent paralogs (70–90% identity):

```bash
crossgene \
  --gene-a HSP90AB1 \
  --gene-b HSP90AA1 \
  --sensitivity 1 \
  --stringency 1
```

### With flanking regions and blacklist masking

Include 5 kb flanking regions and mask repetitive elements before alignment:

```bash
crossgene \
  --gene-a BRCA1 \
  --gene-b BRCA2 \
  --flanking 5000 \
  --blacklist repeats.bed
```

### With BED annotation overlay

Overlay up to 3 BED files on the circular plot as inner annotation rings (e.g., repeat elements, regulatory regions):

```bash
crossgene \
  --gene-a BRCA1 \
  --gene-b BRCA2 \
  --bed repeats.bed \
  --bed regulatory.bed \
  --bed-color darkorchid \
  --bed-color teal
```

Each `--bed` file gets its own ring below the tick track. Colors are auto-assigned from a default palette (darkorchid, teal, sienna) or overridden with `--bed-color`. Region names from BED column 4 are drawn as labels on sufficiently large regions.

### With custom filtering

```bash
crossgene \
  --gene-a BRCA1 \
  --gene-b BRCA2 \
  --min-quality 50 \
  --min-length 40 \
  --max-evalue 0.01 \
  --outdir results/ \
  --output-formats bigwig,tsv,plot \
  --verbose
```

### All parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Gene A | `--gene-a` | required | First gene name (e.g., BRCA1) |
| Gene B | `--gene-b` | required | Second gene name (e.g., BRCA2) |
| Sensitivity | `--sensitivity` | `3` | BLAST sensitivity level (1=most sensitive, 5=strictest) |
| Stringency | `--stringency` | `3` | Post-alignment filtering stringency (1=permissive, 5=strictest) |
| Min quality | `--min-quality` | per stringency | Minimum identity (0–100) to report a hit |
| Min MAPQ | `--min-mapq` | per stringency | Minimum pseudo-MAPQ (derived from bitscore) |
| Min length | `--min-length` | per stringency | Minimum HSP alignment length |
| Min bitscore | `--min-bitscore` | per stringency | Minimum bitscore to report a hit |
| Max E-value | `--max-evalue` | per stringency | Maximum E-value to report a hit |
| Max secondary | `--max-secondary` | `10` | Max HSPs per gene pair |
| Flanking | `--flanking` | `2000` | Flanking region size in bp (upstream + downstream) |
| Blacklist | `--blacklist` | none | BED file of regions to mask with N's before alignment |
| Genome FASTA | `--genome` | `references/homo_sapiens.109.mainChr.fa` | Path to genome FASTA |
| Gene GTF | `--genes-gtf` | `references/homo_sapiens.109.genes.gtf` | Path to gene annotation GTF |
| Annotation GTF | `--annotation-gtf` | `references/homo_sapiens.109.mainChr.gtf` | GTF with sub-gene features (exons, CDS, UTRs) |
| Annotation features | `--annotation-features` | `exon,CDS` | Comma-separated feature types to display |
| Transcript mode | `--transcript-mode` | `canonical` | Transcript selection: `canonical` or `all` |
| Chrom sizes | `--chrom-sizes` | `references/homo_sapiens.109.chrom.sizes` | Chromosome sizes file |
| Output dir | `--outdir` | `.` | Output directory |
| Output formats | `--output-formats` | `bigwig,tsv,plot,bed` | Comma-separated: `bigwig`, `tsv`, `plot`, `bed` |
| BED primary only | `--bed-primary-only` | off | Only include primary HSPs in hit BED export |
| BED overlay | `--bed` | none | BED file for annotation overlay on circular plot (repeatable, max 3) |
| BED color | `--bed-color` | auto | Color for corresponding `--bed` file (default: darkorchid, teal, sienna) |
| Verbose | `--verbose` / `-v` | off | Enable debug logging |

### Sensitivity and stringency presets

The `--sensitivity` and `--stringency` flags provide preset parameter bundles on a 1–5 scale.

**Sensitivity** controls BLAST alignment parameters (word size, scoring, E-value threshold):

| Level | Word size | Reward/Penalty | Gap open/extend | E-value | Use case |
|-------|-----------|----------------|-----------------|---------|----------|
| 1 | 7 | 1 / -1 | 2 / 1 | 10 | High divergence / paralogs (~70–80%) |
| 2 | 7 | 1 / -1 | 2 / 1 | 1 | Moderate divergence (~80–90%) |
| 3 | 11 | 1 / -2 | 1 / 1 | 0.01 | Default — balanced |
| 4 | 13 | 1 / -3 | 1 / 2 | 0.001 | High similarity (~95%+) |
| 5 | 15 | 1 / -4 | 1 / 2 | 0.001 | Near-identical sequences |

**Stringency** controls post-alignment filtering thresholds:

| Level | Min quality | Min MAPQ | Min length | Min bitscore | Max E-value | Use case |
|-------|-------------|----------|------------|--------------|-------------|----------|
| 1 | 0 | 0 | 15 | 0 | inf | Permissive — keep everything |
| 2 | 15 | 0 | 20 | 0 | inf | Mild filtering |
| 3 | 30 | 0 | 25 | 0 | inf | Default — balanced |
| 4 | 50 | 3 | 40 | 30 | 0.01 | Strict |
| 5 | 70 | 5 | 50 | 50 | 0.001 | Very strict — high-confidence only |

Explicitly provided filter flags (e.g., `--min-quality 50`) override the corresponding stringency preset value.

## How It Works

### Pipeline

1. **Gene extraction** — Looks up gene coordinates in the GTF by `gene_name`, extracts the sequence from the genome FASTA using pysam. Sequences are stored in genomic orientation (no reverse-complement). Optional flanking regions are added upstream and downstream.

2. **Blacklist masking** (optional) — Replaces bases in blacklisted regions with N's before alignment.

3. **BLASTN alignment** — Aligns the full query gene sequence against the full target gene sequence using BLASTN. BLAST natively finds all high-scoring pairs (HSPs) without requiring fragmentation.

4. **Parsing** — Parses BLASTN tabular output and translates sequence-local coordinates back to genomic coordinates. Pseudo-MAPQ is derived from bitscore: `min(60, bitscore / max_bitscore * 60)`. Relative strand is determined from both genes' original strands and the BLAST alignment direction.

5. **Filtering** — Applies quality, length, bitscore, E-value, and MAPQ thresholds (from stringency preset or explicit flags).

6. **Scoring** — For each base in the query gene, computes a similarity score (0–100) as the maximum HSP identity across all overlapping hits.

7. **Output** — Writes BigWig tracks, TSV hit tables, BED9 hit files, and a circular visualization.

Steps 3–7 run in both directions (A→B and B→A). The circular plot uses A→B hits only.

### Similarity score

The per-base similarity score (0–100) is computed as:

- For each base position, collect all HSPs that overlap it
- Take the maximum identity among overlapping HSPs
- Scale to 0–100

A score of **100** means perfect match, **0** means no HSP covered that base.

### Visualization

The circular plot (pycirclize) shows:

- **Gene tracks** (outer ring) — Gene A on the right, Gene B on the left, with exon/CDS/UTR annotations drawn as colored rectangles (exons: orange, CDS: forestgreen, 5' UTR: mediumpurple, 3' UTR: goldenrod). Flanking regions shown in lavender if present.
- **Tick track** — Genomic position labels with interval auto-scaled to gene length (1/5/20 kbp).
- **BED annotation rings** (optional, up to 3) — Inner rings below the tick track, each showing regions from a `--bed` file as colored rectangles with labels.
- **Arcs** — Connect matching HSP regions between genes. Color: steelblue (same-sense) or firebrick (antisense). Transparency scales with identity and hit density. Auto-subsampled to 5000 arcs if there are more (top by bitscore).

## Development

### Running tests

```bash
conda activate crossgene
pytest tests/ -v
```

### Project structure

```
crossgene/
├── crossgene/
│   ├── __init__.py
│   ├── cli.py              # Click CLI entry point
│   ├── gene_extractor.py   # GTF lookup + FASTA extraction + feature loading
│   ├── blastn.py            # BLASTN wrapper: alignment, sequence masking
│   ├── parser_blast.py      # BLAST tabular parsing + coordinate translation
│   ├── scores.py            # Per-base scoring (max identity per base)
│   ├── bigwig.py            # BigWig writer (real genomic coordinates)
│   ├── tsv_writer.py        # TSV hit table output
│   ├── bed_writer.py        # BED9 hit export with cross-reference numbering
│   ├── visualize.py         # Circular plot (pycirclize + matplotlib)
│   ├── bed_parser.py        # BED file parsing + region filtering/clipping
│   ├── presets.py           # Sensitivity + stringency preset loader
│   ├── presets.json         # Preset parameter definitions
│   └── models.py            # Data structures (GeneRecord, AlignmentHit, etc.)
├── tests/                   # pytest test suite (one file per module)
├── references/              # Reference data (not tracked)
├── artifacts/               # Design documents and implementation plans
└── pyproject.toml
```

## License

MIT
