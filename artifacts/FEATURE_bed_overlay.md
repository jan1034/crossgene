# Feature: BED File Overlay on Circular Plot

## Problem Statement

Users want to overlay custom genomic annotations (e.g., repeat elements, regulatory regions, variants) from BED files onto the circular comparison plot. Regions should appear as colored blocks on the gene sectors.

**BED format**: `chr1  10223524  10223843  AluSc  2224  -` (6-column BED)

### Constraints
- BED regions must be clipped to each gene's genomic range (only overlapping portions shown)
- Up to 3 BED files supported (each gets its own annotation ring)
- Must not break existing visualization when no BED files are provided
- BED regions may be on either strand or unstranded
- Regions are shown on both Gene A and Gene B sectors wherever they overlap

---

## Current Visualization Layout

```
r=115         Sector labels ("BRCA1 (+)")
r=95-100      Gene body track (exons, CDS, UTRs, flanking)
r=90-95       Tick marks / coordinate labels
r=0-90        Arc space (alignment links between sectors)
```

Two sectors: Gene A (right), Gene B (left), separated by 15° gaps.

---

## Decided Design: Inner Annotation Rings

BED tracks are placed as inner rings below the tick track, above the arc space:

```
r=115         Sector labels
r=95-100      Gene body track (existing)
r=90-95       Tick track (existing)
r=85-90       BED track 1 (e.g., repeats)        ← NEW
r=80-85       BED track 2 (e.g., regulatory)      ← NEW (optional)
r=75-80       BED track 3                          ← NEW (optional)
r=0-75        Arc space
```

**Rationale:**
- Clean visual separation from gene features
- Each BED file gets its own ring → no overlap ambiguity
- Arcs still have plenty of radial space (75px minimum)
- Intuitive reading order: gene structure → annotations → comparisons

---

## Design Decisions

| Question | Decision |
|----------|----------|
| **Placement** | Inner rings (below ticks, above arcs) |
| **Multiple BED files** | Each file gets its own ring, max 3 rings |
| **Coloring** | Single color per BED file (auto-assigned from palette, user can override) |
| **Gene assignment** | Auto-detect: regions shown on both sectors wherever they overlap (by chromosome match) |
| **Labels** | Text labels displayed on the plot (region name from BED column 4) |
| **CLI** | `--bed repeats.bed` (repeatable: `--bed a.bed --bed b.bed --bed c.bed`) |

---

## Color Palette

Auto-assigned colors for up to 3 BED tracks (user can override via `--bed-color`):

| Track | Default Color |
|-------|--------------|
| BED 1 | `darkorchid` |
| BED 2 | `teal` |
| BED 3 | `sienna` |

These are chosen to be visually distinct from existing feature colors (orange, green, purple, gold) and alignment arc colors (steelblue, firebrick).

---

## Data Flow

```
CLI: --bed repeats.bed --bed regulatory.bed
         ↓
    bed_parser.py
    parse each BED file → list[BedRegion]
         ↓
    filter & clip to gene_a range → regions_a
    filter & clip to gene_b range → regions_b
         ↓
    visualize.py
    for each BED file:
      add inner track per sector
      draw rect() blocks for overlapping regions
      draw text labels (region name from col 4)
      add legend entry
```

---

## New Dataclass (`models.py`)

```python
@dataclass
class BedRegion:
    chrom: str
    start: int       # 0-based
    end: int         # exclusive
    name: str        # column 4 (e.g., "AluSc"), default "." if missing
    score: int       # column 5, default 0 if missing
    strand: str      # column 6 ('+', '-', or '.'), default '.' if missing
```

---

## New Module: `bed_parser.py`

Responsibilities:
- Parse standard BED format (tolerant of BED3 through BED6, missing columns get defaults)
- Skip comment lines (`#`) and empty lines
- Validate coordinates (start < end, start >= 0)
- Log warnings for malformed lines (skip, don't crash)

Public functions:

```python
def parse_bed(path: str) -> list[BedRegion]:
    """Parse a BED file into a list of BedRegion objects."""

def filter_and_clip(regions: list[BedRegion], chrom: str, start: int, end: int) -> list[BedRegion]:
    """Return regions overlapping [start, end) on chrom, clipped to boundaries."""
```

---

## Visualization Changes (`visualize.py`)

### New Config Dataclass

```python
@dataclass
class BedTrackConfig:
    label: str                  # Legend label (derived from filename)
    color: str                  # Track color
    regions_a: list[BedRegion]  # Regions overlapping Gene A (already clipped)
    regions_b: list[BedRegion]  # Regions overlapping Gene B (already clipped)
```

### Modified `create_circlize_plot()` Signature

```python
def create_circlize_plot(
    hits_ab: list[AlignmentHit],
    gene_a: GeneRecord,
    gene_b: GeneRecord,
    output_path: str,
    max_arcs: int = 5000,
    bed_tracks: list[BedTrackConfig] | None = None,  # ← NEW
) -> None:
```

### Track Drawing Logic

For each `BedTrackConfig` (index `i` from 0 to 2):
1. Compute radial range: `r_low = 85 - i*5`, `r_high = 90 - i*5`
2. For each sector, add track at `(r_low, r_high)`
3. Draw `track.rect(local_start, local_end, fc=color, ec="none")` for each region
4. Draw `track.text(name, position, size=5)` for region labels (centered on region)
5. Add legend entry: colored patch with `label`

### Label Strategy

- Region names (BED column 4) are drawn as text centered on each region block
- Use `track.text()` at the midpoint of each region
- Font size: 5pt (small to avoid clutter)
- For very small regions where text would overflow, skip the label (threshold: region < 1% of gene length)

---

## CLI Changes (`cli.py`)

```python
@click.option(
    "--bed",
    multiple=True,
    type=click.Path(exists=True),
    help="BED file for annotation overlay on plot (repeatable, max 3).",
)
@click.option(
    "--bed-color",
    multiple=True,
    type=str,
    help="Color for corresponding --bed file (default: auto-assigned).",
)
```

Validation:
- Error if more than 3 `--bed` files provided
- Error if `--bed-color` count exceeds `--bed` count
- Warning if BED file has zero regions overlapping either gene

---

## Testing Strategy

### `test_bed_parser.py`
- Parse valid BED6, BED5, BED4, BED3 files
- Skip comment/empty lines
- Malformed line handling (warning, not crash)
- `filter_and_clip`: overlap logic, boundary clipping, chromosome filtering
- Edge cases: region entirely outside gene, region partially overlapping, region fully contained

### `test_visualize.py` (additions)
- Plot with 1, 2, 3 BED tracks renders without error
- Plot with 0 BED tracks unchanged (backward compat)
- Legend includes BED track entries
- Track radii are correct per BED track index

### `test_cli.py` (additions)
- `--bed` option accepts valid file
- Error on > 3 BED files
- `--bed-color` pairing works

---

## References

- [pycirclize Track API](https://moshi4.github.io/pyCirclize/api-docs/track/)
- [pycirclize Circos API](https://moshi4.github.io/pyCirclize/api-docs/circos/)
- [pycirclize Plot Tips (grouping/rect)](https://moshi4.github.io/pyCirclize/plot_tips/)
