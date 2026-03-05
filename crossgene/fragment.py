"""Rolling-window fragmentation of gene sequences."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

from crossgene.models import BedRegion, GeneRecord

logger = logging.getLogger(__name__)


def _overlaps_blacklist(
    genomic_start: int, genomic_end: int,
    blacklist: list[BedRegion], bl_idx: int,
) -> tuple[bool, int]:
    """Check if a genomic interval overlaps any blacklist region.

    Uses a sweep pointer (bl_idx) for O(n+m) performance with sorted blacklist.
    Advances bl_idx past regions that end before genomic_start.
    Returns (overlaps, updated_bl_idx).
    """
    # Advance past blacklist regions that end at or before our start
    while bl_idx < len(blacklist) and blacklist[bl_idx].end <= genomic_start:
        bl_idx += 1

    # Check remaining regions that could overlap
    i = bl_idx
    while i < len(blacklist) and blacklist[i].start < genomic_end:
        if blacklist[i].end > genomic_start:
            return True, bl_idx
        i += 1

    return False, bl_idx


def generate_fragments(
    gene: GeneRecord, fragment_size: int, step_size: int,
    blacklist: list[BedRegion] | None = None,
) -> Path:
    """Fragment a gene's sense sequence using a rolling window.

    Writes fragments to a temporary FASTA file with headers ``>gene_name:idx``
    where idx is a 0-based sequential index.

    Edge cases:
    - Sequence shorter than fragment_size: single fragment = full sequence.
    - Last fragment shorter than fragment_size: included if its length
      >= fragment_size / 2, otherwise skipped.

    If blacklist is provided, fragments overlapping any blacklisted region
    (in genomic coordinates) are skipped entirely.

    Returns the path to the temporary FASTA file.
    """
    seq = gene.sequence
    seq_len = len(seq)

    if seq_len == 0:
        raise ValueError(f"Gene '{gene.name}' has no sequence; call extract_sequence first")

    # Sort blacklist by start coordinate for sweep algorithm
    if blacklist:
        blacklist = sorted(blacklist, key=lambda r: r.start)

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", prefix=f"frags_{gene.name}_", delete=False
    )

    idx = 0
    skipped = 0
    total = 0
    bl_idx = 0  # sweep pointer for blacklist

    if seq_len <= fragment_size:
        total = 1
        if blacklist:
            g_start, g_end = fragment_index_to_genomic(gene, 0, seq_len, step_size)
            hit, bl_idx = _overlaps_blacklist(g_start, g_end, blacklist, bl_idx)
            if hit:
                skipped = 1
            else:
                tmp.write(f">{gene.name}:{idx}\n{seq}\n")
        else:
            tmp.write(f">{gene.name}:{idx}\n{seq}\n")
    else:
        # For minus-strand genes, genomic coords decrease as sense pos increases,
        # so we need separate sweep pointers depending on strand direction.
        # For minus strand, we reverse the blacklist order since genomic coords decrease.
        if blacklist and gene.strand == "-":
            blacklist = sorted(blacklist, key=lambda r: r.start, reverse=True)

        pos = 0
        frag_idx = 0  # sequential counter for index_to_genomic
        while pos + fragment_size <= seq_len:
            total += 1
            write = True
            if blacklist:
                g_start, g_end = fragment_index_to_genomic(gene, frag_idx, fragment_size, step_size)
                if gene.strand == "-":
                    # For minus strand, check overlap directly (no sweep optimization)
                    write = not any(
                        g_start < r.end and g_end > r.start for r in blacklist
                    )
                else:
                    hit, bl_idx = _overlaps_blacklist(g_start, g_end, blacklist, bl_idx)
                    write = not hit
            if write:
                frag = seq[pos : pos + fragment_size]
                tmp.write(f">{gene.name}:{idx}\n{frag}\n")
                idx += 1
            else:
                skipped += 1
            frag_idx += 1
            pos += step_size

        # Handle trailing partial fragment
        if pos < seq_len:
            remainder = seq[pos:]
            if len(remainder) >= fragment_size / 2:
                total += 1
                write = True
                if blacklist:
                    g_start, g_end = fragment_index_to_genomic(gene, frag_idx, fragment_size, step_size)
                    if gene.strand == "-":
                        write = not any(
                            g_start < r.end and g_end > r.start for r in blacklist
                        )
                    else:
                        hit, bl_idx = _overlaps_blacklist(g_start, g_end, blacklist, bl_idx)
                        write = not hit
                if write:
                    tmp.write(f">{gene.name}:{idx}\n{remainder}\n")
                else:
                    skipped += 1

    if blacklist and skipped > 0:
        logger.info(
            "Direction: skipped %d/%d fragments due to blacklist for %s",
            skipped, total, gene.name,
        )

    tmp.close()
    return Path(tmp.name)


def fragment_index_to_genomic(
    gene: GeneRecord, idx: int, fragment_size: int, step_size: int
) -> tuple[int, int]:
    """Convert a fragment index to genomic coordinates (0-based half-open).

    For plus-strand genes:
        sense_start = idx * step_size
        genomic = [gene.start + sense_start, gene.start + sense_start + frag_len)

    For minus-strand genes (sense sequence is reverse-complemented):
        sense_start = idx * step_size
        genomic = [gene.end - sense_start - frag_len, gene.end - sense_start)

    The returned fragment length is capped at the actual sequence length
    to handle the trailing partial fragment correctly.
    """
    seq_len = len(gene.sequence)
    sense_start = idx * step_size
    frag_len = min(fragment_size, seq_len - sense_start)

    if gene.strand == "-":
        genomic_end = gene.end - sense_start
        genomic_start = genomic_end - frag_len
    else:
        genomic_start = gene.start + sense_start
        genomic_end = genomic_start + frag_len

    return genomic_start, genomic_end
