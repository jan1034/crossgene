# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CrossGene is a Python tool that compares two gene sequences by fragmenting one gene, aligning fragments to another using minimap2, and producing similarity profiles (BigWig), detailed hit tables (TSV), and circular visualizations (pycirclize). It performs bidirectional comparison (A→B and B→A). Supports BED file overlay on circular plots (`--bed`) for annotation of repeats, regulatory regions, etc.

**Status:** Implementation in progress (Step 1 complete). See artifacts/ANALYSIS.md, artifacts/ARCHITECTURE.md, and artifacts/IMPLEMENTATION.md.

## Development Environment

- Conda environment: `crossgene` (Python 3.11)
- Activate: `conda activate crossgene`
- Package installed in editable mode: `pip install -e .`
- CLI entry point: `crossgene`

## Technology Stack

- Python 3.10+ (IDE configured for 3.11)
- CLI: click
- Gene extraction: pysam (FASTA), custom GTF parsing
- Alignment: minimap2 (external, called via subprocess)
- Scoring: numpy
- BigWig output: pyBigWig
- Visualization: pycirclize + matplotlib
- Data structures: dataclasses
- Testing: pytest (planned)

## Planned Module Architecture

Ten modules in `crossgene/`:
- `cli.py` — Click entry point
- `gene_extractor.py` — GTF lookup + FASTA extraction via pysam
- `fragment.py` — Rolling window fragmentation
- `align.py` — minimap2 subprocess wrapper
- `parser.py` — PAF format parsing + coordinate translation
- `scores.py` — Score aggregation (normalized mappability 0-100)
- `bigwig.py` — BigWig writer with real genomic coordinates (IGV-compatible)
- `tsv_writer.py` — 15-column TSV output
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

- Aligner: minimap2 (better for divergent sequences, native secondary alignments)
- Bidirectional: runs A→B and B→A automatically
- Similarity metric: normalized mappability score (0-100, mean of best-hit identity per base)
- Gene lookup: by gene_name field in GTF, first match if multiple
- Strand handling: extract sense sequence, store original strand info
- minimap2 preset: auto-selected based on fragment_size
- Fragment defaults: size=50bp, step=1bp

## Workflow

Per user instructions, follow the phased workflow:
1. **Analysis** → ANALYSIS.md (done)
2. **Architecture** → ARCHITECTURE.md (done)
3. **Implementation plan** → IMPLEMENTATION.md (done)
4. **Implementation** — step by step with user review between steps (Step 1 complete)
