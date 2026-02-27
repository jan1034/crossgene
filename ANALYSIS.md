# Compare Genes — Initial Draft

## Problem Statement

Build a tool that compares the sequences of two genes (Gene A vs Gene B) by:
1. Fragmenting Gene A with a rolling window of configurable size
2. Aligning each fragment to Gene B using an aligner (e.g., BWA)
3. Tracing back fragment alignments to produce similarity profiles and visualizations

### Constraints
- Input: Two gene sequences (FASTA)
- Output: Similarity tracks, coordinate mappings, and visualizations

---

## Core Workflow

```
Gene A (FASTA)          Gene B (FASTA)
     |                       |
  Rolling window          BWA index
  (fragment size, step)      |
     |                       |
  Fragments (FASTA)  ------> BWA alignment
                             |
                        BAM / SAM
                             |
              +--------------+--------------+
              |              |              |
         BigWig track    TSV mapping    Circlize plot
         (similarity     (coord A →     (pycirclize,
          per position)   best match B)  colored arcs)
```

---

## Proposed Outputs

| Output | Format | Description |
|--------|--------|-------------|
| Similarity track | BigWig | Per-base similarity score along Gene A (derived from alignment quality/identity of overlapping fragments) |
| Coordinate mapping | TSV | Each fragment's position in A → best hit position in B, with alignment score, identity %, strand |
| Visualization | PDF/SVG | Pycirclize plot: Gene A on one side, Gene B on the other, arcs connecting matched regions colored by alignment quality |

---

## Key Parameters (initial list)

| Parameter | Description | Example default |
|-----------|-------------|-----------------|
| `fragment_size` | Length of each rolling-window fragment | 100 bp |
| `step_size` | Step between consecutive fragments | 1 bp (or fraction of fragment_size) |
| `min_quality` | Minimum MAPQ or identity to report a match | 30 |
| `aligner` | Alignment tool | `bwa mem` |

---

## Open Design Questions

### 1. Fragment generation
- **Step size**: Should default be 1 bp (maximum resolution, slow) or e.g. `fragment_size / 2` for overlap?
- Should we support variable fragment sizes in one run?

### 2. Alignment
- **BWA mem vs minimap2?** minimap2 may be better for longer fragments or divergent sequences.
- How to handle multiple alignments per fragment? (multimappers)
  - Report best only? Top N? All above threshold?
- Should we consider reverse complement matches?

### 3. Similarity score computation
- How to aggregate overlapping fragment scores into a per-base similarity value?
  - Options: mean alignment identity, max identity, weighted by MAPQ
- What metric: % identity, alignment score, e-value equivalent?

### 4. BigWig generation
- Need chromosome sizes — since these are individual genes, we'd treat each gene as a "chromosome" with length = gene length.
- Tool: `deeptools` or `pyBigWig` for writing BigWig files.

### 5. TSV output
- Columns to include? Proposal:
  ```
  frag_start_A  frag_end_A  hit_start_B  hit_end_B  strand  identity  mapq  cigar
  ```

### 6. Visualization
- **Pycirclize**: Gene A as right semicircle, Gene B as left semicircle, Bezier arcs connecting matched regions.
- Color scale: similarity (e.g., red=high, blue=low) or discrete bins?
- Should we add gene feature annotations (exons, domains) if provided?
- Filter: only show top N arcs or arcs above a quality threshold to avoid clutter?

### 7. Scope
- Gene-level only, or should this scale to whole chromosomes / genomes?
- Single pair comparison only, or batch mode?

---

## Technology Stack (tentative)

| Component | Tool/Library |
|-----------|-------------|
| Language | Python (CLI entry point) |
| Fragmentation | Custom (BioPython or plain Python) |
| Alignment | BWA mem (subprocess) or minimap2 |
| BAM parsing | pysam |
| BigWig writing | pyBigWig |
| TSV | pandas or plain csv |
| Visualization | pycirclize |
| CLI | argparse or click |

---

## Next Steps

1. Discuss and resolve open design questions above
2. Produce `ARCHITECTURE.md` with module boundaries and data flow
3. Write `IMPLEMENTATION.md` with step-by-step plan
4. Implement incrementally with review checkpoints
