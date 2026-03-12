# Debugging Plan: Blacklist & Sensitivity Issues

## Issue 1: `--blacklist` Not Recognized

### Status: Cannot reproduce from code

The CLI definition in `cli.py:130` is correct:
```python
@click.option("--blacklist", default=None, type=click.Path(exists=True), ...)
```

The test script `process_2x.sh` ran successfully and produced output in all directories including blacklisted runs (`no-sines-*`). The blacklist correctly reduces hit extent (e.g., strict mode: 170bp hit without blacklist vs 141bp with blacklist on the primary CDS alignment).

### Root Cause Hypothesis

- **Stale installation:** The installed CLI entry point may not reflect the latest code. If the package was installed with `pip install -e .` but the conda environment was rebuilt or the entry point was cached, the old version without `--blacklist` could be running.

### Action Items

1. Verify installed version: `which crossgene` and `pip show crossgene`
2. Reinstall: `pip install -e .` in the `compare_genes` conda environment
3. Confirm: `crossgene --help | grep blacklist`

---

## Issue 2: Hit Count Independent of Sensitivity Mode

### Status: Confirmed — all non-strict modes produce exactly 10 hits

| Mode | Hits (A→B) | Notes |
|------|-----------|-------|
| strict | 1 | max_secondary=1, min_quality=70 |
| default | 10 | capped at max_secondary |
| moderate | 10 | capped at max_secondary |
| sensitive | 10 | capped at max_secondary |
| divergent | 10 | capped at max_secondary |

### Root Cause: `max_hsps` Hard Cap

In `blastn.py:83`:
```python
"-max_hsps", str(params.max_secondary),  # default=10
```

BLAST's `-max_hsps` parameter limits the number of HSPs reported **regardless** of how many exist. All sensitivity modes find >=10 HSPs for Actg1/Actg2 (highly similar actin paralogs), so the output is capped at 10 for all of them.

The sensitivity parameters DO change BLAST behavior (different word_size, reward/penalty produce different hit coordinates and identities), but the **count** appears identical because of the cap.

### Secondary Problem: Sensitivity Presets Are Too Similar

Comparing divergent vs sensitive mode — they use identical parameters except evalue (10 vs 1):
```
divergent: word_size=7, reward=1, penalty=-1, gapopen=2, gapextend=1, evalue=10
sensitive: word_size=7, reward=1, penalty=-1, gapopen=2, gapextend=1, evalue=1
```

This means for most gene comparisons, divergent and sensitive produce nearly identical results.

---

## Proposed Fix: Replace Flags with `--sensitivity 1-10`

### Design

Replace `--moderate`, `--sensitive`, `--divergent`, `--strict` with a single `--sensitivity` integer parameter (1-10) that provides a smooth gradient of BLAST tuning.

### Parameter Mapping

| Parameter | 1 (strictest) | 3 | 5 (default) | 7 | 10 (most sensitive) |
|-----------|--------------|---|-------------|---|---------------------|
| word_size | 28 | 18 | 11 | 7 | 7 |
| reward | 1 | 1 | 1 | 1 | 1 |
| penalty | -5 | -3 | -2 | -1 | -1 |
| gapopen | 2 | 2 | 1 | 2 | 0 |
| gapextend | 2 | 1 | 1 | 1 | 2 |
| evalue | 1e-10 | 1e-5 | 0.01 | 1 | 100 |
| max_hsps | 5 | 10 | 25 | 50 | 100 |
| min_quality | 70 | 50 | 30 | 20 | 0 |
| min_length | 50 | 30 | 25 | 15 | 10 |

Key insight: **`max_hsps` must scale with sensitivity** to actually affect hit counts. Currently it's fixed at 10 regardless of mode.

### Implementation Steps

#### Step 1: Define parameter interpolation function

Create a `_sensitivity_params(level: int) -> dict` function in `blastn.py` that interpolates BLAST parameters across the 1-10 range. Use discrete breakpoints (not linear interpolation) since BLAST parameters interact non-linearly.

```python
# Breakpoints: (level, word_size, penalty, gapopen, gapextend, evalue, max_hsps)
_PRESETS = [
    (1,  28, -5, 2, 2, 1e-10, 5),
    (2,  22, -4, 2, 2, 1e-7,  5),
    (3,  18, -3, 2, 1, 1e-5,  10),
    (4,  15, -3, 1, 1, 1e-3,  15),
    (5,  11, -2, 1, 1, 0.01,  25),
    (6,  11, -2, 1, 1, 0.1,   35),
    (7,   7, -1, 2, 1, 1,     50),
    (8,   7, -1, 2, 1, 5,     75),
    (9,   7, -1, 0, 2, 10,    100),
    (10,  7, -1, 0, 2, 100,   100),
]
```

#### Step 2: Update `BlastParams` dataclass

```python
@dataclass
class BlastParams:
    sensitivity: int = 5    # 1-10
    max_hsps: int = 25      # derived from sensitivity
    # Remove: moderate, sensitive, divergent booleans
```

#### Step 3: Update CLI

```python
# Remove these:
# --moderate, --sensitive, --divergent, --strict

# Add:
@click.option("--sensitivity", type=click.IntRange(1, 10), default=5,
              help="BLAST sensitivity 1-10 (1=strict/fast, 10=divergent/slow)")
```

#### Step 4: Update `_build_blastn_command()`

Replace the if/elif chain with a lookup into the preset table.

#### Step 5: Update filter defaults

Make `min_quality` and `min_length` also sensitivity-aware (lower thresholds at higher sensitivity), unless the user explicitly overrides them.

#### Step 6: Update tests

Update any tests that use the old flag names.

### Migration Notes

- `--strict` ≈ `--sensitivity 1`
- default ≈ `--sensitivity 5`
- `--moderate` ≈ `--sensitivity 6`
- `--sensitive` ≈ `--sensitivity 7`
- `--divergent` ≈ `--sensitivity 9`

### Validation Plan

After implementation, re-run `process_2x.sh` (adapted for new syntax) and verify:
1. Hit counts increase monotonically with sensitivity level
2. Actg1 vs Actg2 (similar paralogs): sensitivity 1 → few hits, sensitivity 10 → many hits
3. Lmna vs Hprt (unrelated): even sensitivity 10 produces few/no hits
4. Blacklist reduces hits consistently across all levels
5. BigWig coverage increases with sensitivity
