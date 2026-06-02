# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CrossGene is a Python tool that compares two gene sequences using direct BLASTN alignment (gene-vs-gene). It produces similarity profiles (BigWig), detailed hit tables (TSV), BED files, and circular visualizations (pycirclize). It performs bidirectional comparison (A→B and B→A). Supports BED file overlay on circular plots (`--bed`) for annotation of repeats, regulatory regions, etc.

## Development Environment

- Conda environment: `crossgene` (Python 3.11)
- Activate: `conda activate crossgene`
- Package installed in editable mode: `pip install -e .`
- CLI entry point: `crossgene`

## Technology Stack

- Python 3.10+ (IDE configured for 3.11)
- CLI: click
- Gene extraction: pysam (FASTA), custom GTF parsing
- Alignment: BLASTN (BLAST+, called via subprocess)
- Scoring: numpy
- BigWig output: pyBigWig
- Visualization: pycirclize + matplotlib
- Data structures: dataclasses
- Testing: pytest

## Module Architecture

Modules in `crossgene/`:
- `cli.py` — Click entry point
- `gene_extractor.py` — GTF lookup + FASTA extraction via pysam
- `blastn.py` — BLASTN wrapper: gene-vs-gene alignment, sequence masking
- `parser_blast.py` — BLAST tabular output parser + coordinate translation
- `scores.py` — HSP-based per-base scoring (best identity per base, 0-100)
- `bigwig.py` — BigWig writer with real genomic coordinates (IGV-compatible)
- `tsv_writer.py` — TSV output with evalue/bitscore columns
- `bed_writer.py` — BED9 output for alignment hits
- `visualize.py` — Circlize circular plot
- `models.py` — GeneRecord, GeneFeature, AlignmentHit, BedRegion dataclasses
- `bed_parser.py` — BED file parsing + region filtering/clipping

## Reference Data

Located in `references/` (gitignored, ~2.1 GB):
- `homo_sapiens.109.mainChr.fa` — Human genome FASTA (Ensembl 109, main chromosomes)
- `homo_sapiens.109.mainChr.fa.fai` — FASTA index
- `homo_sapiens.109.genes.gtf` — Gene annotations
- `homo_sapiens.109.chrom.sizes` — Chromosome sizes for BigWig generation

## Key Design Decisions

- Aligner: BLASTN gene-vs-gene (no fragmentation, BLAST finds HSPs natively)
- Bidirectional: runs A→B and B→A automatically
- Similarity metric: per-base best HSP identity (0-100, max wins for overlapping HSPs)
- Gene lookup: by gene_name field in GTF, first match if multiple
- Strand handling: extract sense sequence, store original strand info
- Sensitivity modes: `--moderate`, `--sensitive`, `--divergent` (word_size/scoring presets)
- Blacklist: hard-mask regions with N's before BLAST alignment
- Filters: `--min-quality`, `--min-length`, `--min-bitscore`, `--max-evalue`

## Pipeline Flow

```
Gene A sequence (full)
    │
    ▼
BLASTN: gene A (query) → gene B (subject)
    │
    ▼
Parse HSPs → list[AlignmentHit]
    │
    ├─► Filter by min_quality, min_length, max_evalue, min_bitscore
    │
    ├─► compute_scores() → BigWig (per-base best identity)
    ├─► write_tsv()
    ├─► write_hit_beds()
    └─► circlize plot

(repeat B→A direction)
```
