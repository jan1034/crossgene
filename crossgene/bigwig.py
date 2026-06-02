"""BigWig writer for per-base mappability scores."""

from __future__ import annotations

import numpy as np
import pyBigWig


def read_chrom_sizes(chrom_sizes_path: str) -> dict[str, int]:
    """Read a chrom.sizes file into a dict of {chrom: length}.

    Handles optional header lines (starting with non-chr or containing
    alphabetic characters in the second column).
    """
    sizes: dict[str, int] = {}
    with open(chrom_sizes_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError:
                # Skip header or malformed lines
                continue
    return sizes


def write_bigwig(
    scores: np.ndarray,
    chrom: str,
    start: int,
    chrom_sizes: dict[str, int],
    output_path: str,
) -> None:
    """Write per-base scores as a BigWig file with genomic coordinates.

    Each base position ``start + i`` gets ``scores[i]``. Uses batch
    ``addEntries`` for efficiency.

    Args:
        scores: Per-base score array (0-100), length = gene.end - gene.start.
        chrom: Chromosome name (e.g. "chr17").
        start: Genomic start position (0-based).
        chrom_sizes: Dict of {chrom: length}.
        output_path: Output BigWig file path.
    """

    if chrom not in chrom_sizes:
        raise ValueError(
            f"Chromosome '{chrom}' not found in chrom.sizes file. "
            f"Available: {', '.join(sorted(chrom_sizes.keys()))}"
        )

    bw = pyBigWig.open(output_path, "w")
    try:
        # Add all chromosome headers (required by pyBigWig)
        headers = list(chrom_sizes.items())
        bw.addHeader(headers)

        # Write scores as fixed-step entries (one value per base)
        n = len(scores)
        if n == 0:
            return

        chroms = np.array([chrom] * n)
        starts = np.arange(start, start + n, dtype=np.int64)
        ends = starts + 1
        values = scores.astype(np.float32)

        bw.addEntries(chroms, starts, ends=ends, values=values)
    finally:
        bw.close()
