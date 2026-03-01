# BED Overlay — Implementation Plan

Based on decisions in [FEATURE_bed_overlay.md](./FEATURE_bed_overlay.md).

---

## Steps

### Step 1 — BedRegion Dataclass (`models.py`)

Add the `BedRegion` dataclass to the existing models module.

**Changes to `crossgene/models.py`:**

```python
@dataclass
class BedRegion:
    """A region from a BED file."""
    chrom: str
    start: int       # 0-based
    end: int         # exclusive
    name: str        # column 4, default "."
    score: int       # column 5, default 0
    strand: str      # column 6, default "."
```

**Checkpoint:** Import and instantiate `BedRegion` in a Python shell.

---

### Step 2 — BED Parser (`bed_parser.py`)

Create a new module for parsing BED files and filtering/clipping regions to gene boundaries.

**File:** `crossgene/bed_parser.py`

**Functions:**

```python
def parse_bed(path: str) -> list[BedRegion]:
    """Parse a BED file (BED3 through BED6).

    - Skip lines starting with '#', 'track', 'browser', or empty lines
    - Missing columns get defaults: name=".", score=0, strand="."
    - Log warning and skip malformed lines (start >= end, non-numeric coords)
    """

def filter_and_clip(
    regions: list[BedRegion], chrom: str, start: int, end: int
) -> list[BedRegion]:
    """Return regions overlapping [start, end) on chrom, clipped to boundaries.

    - Filter: region.chrom == chrom AND region.start < end AND region.end > start
    - Clip: region.start = max(region.start, start), region.end = min(region.end, end)
    - Returns new BedRegion objects (does not mutate input)
    """
```

**Tests (`tests/test_bed_parser.py`):**
- Parse valid BED6 file → correct BedRegion list
- Parse BED3 file → defaults for name, score, strand
- Skip comment/track/browser/empty lines
- Malformed lines: warning logged, line skipped, no crash
- `filter_and_clip`: region fully inside → unchanged
- `filter_and_clip`: region partially overlapping → clipped
- `filter_and_clip`: region fully outside → excluded
- `filter_and_clip`: wrong chromosome → excluded
- `filter_and_clip`: edge case start == end after clipping → excluded

**Checkpoint:** Parse a sample BED6 file, filter to a known genomic interval, verify output.

---

### Step 3 — BedTrackConfig and Visualization Changes (`visualize.py`)

Add BED track rendering to the circular plot.

**New dataclass** (in `visualize.py` or `models.py`, keep local to visualize since it's visualization-specific):

```python
@dataclass
class BedTrackConfig:
    label: str                  # legend label (derived from filename stem)
    color: str                  # single color for all regions in this track
    regions_a: list[BedRegion]  # regions overlapping Gene A (pre-clipped)
    regions_b: list[BedRegion]  # regions overlapping Gene B (pre-clipped)
```

**Changes to `create_circlize_plot()`:**

1. Add optional parameter: `bed_tracks: list[BedTrackConfig] | None = None`
2. After drawing the tick track, for each BED track (index `i`, max 3):
   - Compute radial range: `r_high = 90 - i * 5`, `r_low = r_high - 5`
   - For each sector, add track at `(r_low, r_high)`
   - Draw `track.rect(local_start, local_end, fc=color, ec="none")` per region
   - Draw text label (region name) centered on region via `track.text()` at midpoint, size=5
   - Skip label if region is < 1% of gene length (too small for readable text)
3. Add legend entry per BED track (colored patch with `label`)

**Default color palette:**

```python
BED_COLORS = ["darkorchid", "teal", "sienna"]
```

**Tests (`tests/test_visualize.py` additions):**
- Plot with 1 BED track → PDF created, no errors
- Plot with 3 BED tracks → PDF created, no errors
- Plot with 0 BED tracks (None) → unchanged behavior
- Legend includes BED track entries when present

**Checkpoint:** Generate a plot with synthetic BED regions on both genes, visually inspect.

---

### Step 4 — CLI Integration (`cli.py`)

Wire BED file parsing into the pipeline and pass to visualization.

**New CLI options:**

```python
@click.option(
    "--bed", "bed_files",
    multiple=True,
    type=click.Path(exists=True),
    help="BED file for annotation overlay on circular plot (repeatable, max 3).",
)
@click.option(
    "--bed-color", "bed_colors",
    multiple=True,
    type=str,
    help="Color for corresponding --bed file (default: auto-assigned).",
)
```

**Validation (in `main()`):**
- Error if `len(bed_files) > 3`
- Error if `len(bed_colors) > len(bed_files)`
- Only process BED files when `"plot"` is in output formats (warn otherwise)

**Pipeline changes (before visualization call):**

```python
bed_tracks = []
for i, bed_path in enumerate(bed_files):
    regions = parse_bed(bed_path)
    color = bed_colors[i] if i < len(bed_colors) else BED_COLORS[i]
    label = Path(bed_path).stem  # filename without extension as legend label

    regions_a = filter_and_clip(regions, rec_a.chrom, rec_a.start, rec_a.end)
    regions_b = filter_and_clip(regions, rec_b.chrom, rec_b.start, rec_b.end)

    if not regions_a and not regions_b:
        logger.warning("BED file %s: no regions overlap either gene", bed_path)

    bed_tracks.append(BedTrackConfig(
        label=label, color=color,
        regions_a=regions_a, regions_b=regions_b,
    ))
```

Pass `bed_tracks` (or `None` if empty) to `create_circlize_plot()`.

**Tests (`tests/test_cli.py` additions):**
- `--bed` accepts a valid BED file path
- Error on `--bed` repeated > 3 times
- `--bed-color` pairs with `--bed` correctly
- Warning when BED has no overlapping regions

**Checkpoint:** Full end-to-end run: `crossgene --gene-a X --gene-b Y --bed repeats.bed` produces a plot with annotation rings.

---

### Step 5 — Polish and Edge Cases

- Verify backward compatibility: all existing tests still pass without `--bed`
- Test with real BED files (e.g., RepeatMasker output) on a real gene pair
- Visually check label readability at different gene sizes
- Check behavior with empty BED file (0 regions) → track drawn but empty, no crash
- Check behavior when BED regions only overlap one gene → track appears on one sector only
- Update `CLAUDE.md` project description to mention BED overlay capability

**Checkpoint:** All tests pass. Visual inspection of real-data plot with BED overlay looks correct.

---

## Summary

| Step | Module | Type | Dependencies |
|------|--------|------|-------------|
| 1 | `models.py` | Edit | None |
| 2 | `bed_parser.py` + tests | New file | Step 1 |
| 3 | `visualize.py` + tests | Edit | Step 1 |
| 4 | `cli.py` + tests | Edit | Steps 2, 3 |
| 5 | Polish | All | Steps 1–4 |
