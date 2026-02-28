# Compare Genes

A Python CLI tool for comparing two gene sequences by fragment alignment. It fragments one gene with a rolling window, aligns the fragments to another gene using [minimap2](https://github.com/lh3/minimap2) or [BLASTN](https://blast.ncbi.nlm.nih.gov/), and produces similarity profiles, detailed hit tables, and circular visualizations. The comparison runs bidirectionally (A→B and B→A).

## Outputs

For a comparison of e.g. BRCA1 vs BRCA2, the tool produces:

| File | Format | Description |
|------|--------|-------------|
| `BRCA1_vs_BRCA2.mappability.bw` | BigWig | Per-base similarity score (0–100) of BRCA1 relative to BRCA2 |
| `BRCA2_vs_BRCA1.mappability.bw` | BigWig | Per-base similarity score (0–100) of BRCA2 relative to BRCA1 |
| `BRCA1_vs_BRCA2.hits.tsv` | TSV | All individual alignment hits from A→B |
| `BRCA2_vs_BRCA1.hits.tsv` | TSV | All individual alignment hits from B→A |
| `BRCA1_vs_BRCA2.circlize.pdf` | PDF | Circular plot with arcs connecting matching regions |

BigWig files use real genomic coordinates and can be loaded directly into IGV.

## Installation

### Prerequisites

- Python 3.10+
- [minimap2](https://github.com/lh3/minimap2) on PATH (default aligner)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) on PATH (optional, for `--aligner blastn`)
- Reference data (genome FASTA, GTF annotations, chromosome sizes)

### Install

```bash
conda create -n compare_genes python=3.11
conda activate compare_genes
conda install -c bioconda minimap2 blast

git clone <repository-url>
cd compare_genes
pip install -e .
```

Verify the installation:

```bash
compare-genes --help
```

## Reference Data

The tool requires a genome FASTA, gene annotation GTF, and chromosome sizes file. By default it looks for these in a `references/` directory:

```
references/
├── homo_sapiens.109.mainChr.fa       # Genome FASTA
├── homo_sapiens.109.mainChr.fa.fai   # FASTA index
├── homo_sapiens.109.genes.gtf        # Gene annotations (Ensembl 109)
└── homo_sapiens.109.chrom.sizes      # Chromosome sizes
```

These paths are configurable via CLI flags (see below).

## Usage

### Basic usage

```bash
compare-genes --gene-a BRCA1 --gene-b BRCA2
```

### Using BLASTN (for divergent paralogs)

BLASTN provides better sensitivity for divergent paralogs (70–90% identity) where minimap2 may miss hits:

```bash
compare-genes \
  --gene-a HSP90AB1 \
  --gene-b HSP90AA1 \
  --aligner blastn \
  --divergent
```

### With custom parameters

```bash
compare-genes \
  --gene-a BRCA1 \
  --gene-b BRCA2 \
  --fragment-size 100 \
  --step-size 10 \
  --min-quality 50 \
  --outdir results/ \
  --output-formats bigwig,tsv,plot \
  --verbose
```

### All parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Gene A | `--gene-a` | required | First gene name (e.g., BRCA1) |
| Gene B | `--gene-b` | required | Second gene name (e.g., BRCA2) |
| Fragment size | `--fragment-size` | 50 | Rolling window fragment length in bp |
| Step size | `--step-size` | 1 | Step between fragments in bp |
| Min quality | `--min-quality` | 30 | Minimum identity (0–100) to report a hit |
| Max secondary | `--max-secondary` | 10 | Max secondary alignments per fragment |
| Aligner | `--aligner` | `minimap2` | Alignment engine: `minimap2` or `blastn` |
| Genome FASTA | `--genome` | `references/homo_sapiens.109.mainChr.fa` | Path to genome FASTA |
| GTF | `--gtf` | `references/homo_sapiens.109.genes.gtf` | Path to gene annotation GTF |
| Chrom sizes | `--chrom-sizes` | `references/homo_sapiens.109.chrom.sizes` | Chromosome sizes file |
| Output dir | `--outdir` | `.` | Output directory |
| Output formats | `--output-formats` | `bigwig,tsv,plot` | Comma-separated: `bigwig`, `tsv`, `plot` |
| minimap2 preset | `--minimap2-preset` | `auto` | Override minimap2 preset (auto selects `-x sr` for fragments <=200 bp) |
| Sensitive mode | `--sensitive` | off | Tune aligner for moderate divergence |
| Divergent mode | `--divergent` | off | Tune aligner for high divergence / paralogs (mutually exclusive with `--sensitive`) |
| Verbose | `--verbose` / `-v` | off | Enable debug logging |

## How It Works

### Pipeline

1. **Gene extraction** — Looks up gene coordinates in the GTF by `gene_name`, extracts the sequence from the genome FASTA using pysam. Minus-strand genes are reverse-complemented to produce the sense sequence.

2. **Fragmentation** — Creates overlapping fragments using a rolling window (`--fragment-size` and `--step-size`).

3. **Alignment** — Aligns fragments to the target gene sequence using minimap2 (PAF output) or BLASTN (tabular output).

4. **Parsing** — Parses aligner output and translates fragment-local and target-local coordinates back to genomic coordinates. For BLASTN, MAPQ is derived from bitscore: `min(60, bitscore / max_bitscore * 60)`.

5. **Scoring** — For each base in the query gene, computes a mappability score (0–100) by averaging the best alignment identity across all overlapping fragments.

6. **Output** — Writes BigWig tracks, TSV hit tables, and a circular visualization.

Steps 2–6 run in both directions (A→B and B→A).

### Choosing an aligner

| Aligner | Best for | Notes |
|---------|----------|-------|
| `minimap2` (default) | Highly similar genes, fast runtime | Struggles with divergent paralogs (70–90% identity) |
| `blastn` | Divergent paralogs, homology detection | Better sensitivity at lower identity; slower due to `makeblastdb` step |

Both aligners support `--sensitive` and `--divergent` flags, which adjust internal parameters (seed size, scoring, e-value thresholds) for increased sensitivity. `--minimap2-preset` is ignored when using `--aligner blastn`.

### Similarity score

The per-base mappability score (0–100) is computed as:

- For each base position, collect all fragments that overlap it
- For each fragment, take the best (highest identity) alignment
- Average these best identities across overlapping fragments
- Scale to 0–100

A score of **100** means perfect match, **0** means no fragment aligned at all.

### Visualization

The circular plot (pycirclize) shows:
- Gene A on the right, Gene B on the left, with exon/CDS annotations
- Arcs connecting matching regions (A→B direction)
- Arc color: blue = same-sense match, red = antisense match
- Arc transparency: proportional to mapping quality
- Auto-subsampled to 5000 arcs if there are more

## Development

### Running tests

```bash
conda activate compare_genes
pytest tests/ -v
```

### Project structure

```
compare_genes/
├── compare_genes/
│   ├── __init__.py
│   ├── cli.py              # Click CLI entry point
│   ├── gene_extractor.py   # GTF lookup + FASTA extraction
│   ├── fragment.py          # Rolling window fragmentation
│   ├── align.py             # minimap2 wrapper
│   ├── blastn.py            # BLASTN wrapper
│   ├── parser.py            # PAF parsing + coordinate translation
│   ├── parser_blast.py      # BLAST tabular parsing + coordinate translation
│   ├── scores.py            # Per-base score aggregation
│   ├── bigwig.py            # BigWig writer
│   ├── tsv_writer.py        # TSV output
│   ├── visualize.py         # Circular plot (pycirclize)
│   └── models.py            # Data structures
├── tests/
├── references/              # Reference data (not tracked)
├── pyproject.toml
├── artifacts/
│   ├── ANALYSIS.md
│   ├── ARCHITECTURE.md
│   └── IMPLEMENTATION.md
```

## License

MIT
