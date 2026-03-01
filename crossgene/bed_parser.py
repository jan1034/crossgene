"""BED file parsing and region filtering/clipping."""

from __future__ import annotations

import logging

from crossgene.models import BedRegion

logger = logging.getLogger(__name__)


def parse_bed(path: str) -> list[BedRegion]:
    """Parse a BED file (BED3 through BED6).

    - Skip lines starting with '#', 'track', 'browser', or empty lines
    - Missing columns get defaults: name=".", score=0, strand="."
    - Log warning and skip malformed lines (start >= end, non-numeric coords)
    """
    regions: list[BedRegion] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                logger.warning("Line %d: need at least 3 columns, got %d — skipping", lineno, len(fields))
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                logger.warning("Line %d: non-numeric coordinates — skipping", lineno)
                continue
            if start < 0:
                logger.warning("Line %d: negative start coordinate — skipping", lineno)
                continue
            if start >= end:
                logger.warning("Line %d: start >= end (%d >= %d) — skipping", lineno, start, end)
                continue
            name = fields[3] if len(fields) > 3 else "."
            try:
                score = int(fields[4]) if len(fields) > 4 else 0
            except ValueError:
                score = 0
            strand = fields[5] if len(fields) > 5 else "."
            regions.append(BedRegion(chrom=chrom, start=start, end=end, name=name, score=score, strand=strand))
    return regions


def filter_and_clip(
    regions: list[BedRegion], chrom: str, start: int, end: int
) -> list[BedRegion]:
    """Return regions overlapping [start, end) on chrom, clipped to boundaries.

    - Filter: region.chrom == chrom AND region.start < end AND region.end > start
    - Clip: region.start = max(region.start, start), region.end = min(region.end, end)
    - Returns new BedRegion objects (does not mutate input)
    - Excludes regions where clipped start == clipped end
    """
    result: list[BedRegion] = []
    for r in regions:
        if r.chrom != chrom or r.start >= end or r.end <= start:
            continue
        clipped_start = max(r.start, start)
        clipped_end = min(r.end, end)
        if clipped_start >= clipped_end:
            continue
        result.append(BedRegion(
            chrom=r.chrom, start=clipped_start, end=clipped_end,
            name=r.name, score=r.score, strand=r.strand,
        ))
    return result
