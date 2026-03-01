"""BLASTN tabular output parser with coordinate translation to genomic positions."""

from __future__ import annotations

import logging
from pathlib import Path

from crossgene.fragment import fragment_index_to_genomic
from crossgene.models import AlignmentHit, GeneRecord
from crossgene.parser import _target_local_to_genomic

logger = logging.getLogger(__name__)

# Column indices for -outfmt "6 qseqid sseqid pident length mismatch gapopen
#                               qstart qend sstart send evalue bitscore qlen slen sstrand"
_COL_QSEQID = 0
_COL_PIDENT = 2
_COL_LENGTH = 3
_COL_SSTART = 8
_COL_SEND = 9
_COL_BITSCORE = 11
_COL_SSTRAND = 14
_NUM_COLUMNS = 15


def parse_blast_tabular(
    blast_tsv: Path,
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    fragment_size: int,
    step_size: int,
    min_quality: int,
    direction: str,
    min_mapq: int = 0,
) -> list[AlignmentHit]:
    """Parse BLASTN tabular output into AlignmentHit objects.

    Uses a two-pass approach:
    1. First pass: collect all bitscores to find max (for MAPQ normalization).
    2. Second pass: build AlignmentHit list with normalized mapq.

    Coordinate translation:
    - Query: fragment index → genomic via fragment_index_to_genomic()
    - Target: BLAST sstart/send (1-based inclusive) → 0-based half-open → genomic

    Args:
        blast_tsv: Path to BLASTN tabular output file.
        query_gene: The fragmented (query) gene with genomic coords.
        target_gene: The alignment reference gene with genomic coords.
        fragment_size: Fragment size used for fragmentation.
        step_size: Step size used for fragmentation.
        min_quality: Minimum identity threshold (0-100).
        direction: Label string, e.g. "A→B".

    Returns:
        List of AlignmentHit objects passing the quality filter.
    """
    # Read all lines
    lines: list[str] = []
    with open(blast_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            lines.append(line)

    if not lines:
        logger.debug("Empty BLAST output: %s", blast_tsv)
        return []

    # First pass: find max bitscore for MAPQ normalization
    max_bitscore = 0.0
    for line in lines:
        fields = line.split("\t")
        if len(fields) < _NUM_COLUMNS:
            continue
        bitscore = float(fields[_COL_BITSCORE])
        if bitscore > max_bitscore:
            max_bitscore = bitscore

    if max_bitscore == 0:
        max_bitscore = 1.0  # avoid division by zero

    # Second pass: build hits
    # Track which query names we've seen to assign primary/secondary
    seen_queries: set[str] = set()
    hits: list[AlignmentHit] = []

    for line in lines:
        fields = line.split("\t")
        if len(fields) < _NUM_COLUMNS:
            logger.warning("Skipping malformed BLAST line: %s", line[:80])
            continue

        query_name = fields[_COL_QSEQID]
        pident = float(fields[_COL_PIDENT])  # 0-100 from BLAST
        identity = pident / 100.0

        # Filter by quality
        if pident < min_quality:
            continue

        # BLAST sstart/send: 1-based inclusive
        sstart = int(fields[_COL_SSTART])
        send = int(fields[_COL_SEND])
        bitscore = float(fields[_COL_BITSCORE])
        sstrand = fields[_COL_SSTRAND]

        # Convert strand
        strand = "+" if sstrand == "plus" else "-"

        # Convert sstart/send to 0-based half-open local coordinates
        # BLAST reports sstart > send for minus strand hits
        if sstart <= send:
            local_start = sstart - 1  # 1-based → 0-based
            local_end = send          # inclusive → exclusive
        else:
            local_start = send - 1
            local_end = sstart

        # Derive MAPQ from bitscore
        mapq = min(60, int(bitscore / max_bitscore * 60))

        # Filter by MAPQ
        if mapq < min_mapq:
            continue

        # Primary/secondary: first HSP per query is primary
        is_primary = query_name not in seen_queries
        seen_queries.add(query_name)

        # Extract fragment index from query name (format: "gene_name:idx")
        try:
            frag_idx = int(query_name.rsplit(":", 1)[1])
        except (IndexError, ValueError):
            logger.warning("Cannot parse fragment index from '%s'", query_name)
            continue

        # Translate coordinates to genomic
        query_genomic_start, query_genomic_end = fragment_index_to_genomic(
            query_gene, frag_idx, fragment_size, step_size
        )
        target_genomic_start, target_genomic_end = _target_local_to_genomic(
            target_gene, local_start, local_end
        )

        hits.append(
            AlignmentHit(
                query_chrom=query_gene.chrom,
                query_start=query_genomic_start,
                query_end=query_genomic_end,
                target_chrom=target_gene.chrom,
                target_start=target_genomic_start,
                target_end=target_genomic_end,
                strand=strand,
                identity=identity,
                mapq=mapq,
                cigar="",  # BLAST tabular doesn't include CIGAR
                alignment_score=int(bitscore),
                is_primary=is_primary,
                query_gene=query_gene.name,
                target_gene=target_gene.name,
                direction=direction,
            )
        )

    logger.debug(
        "Parsed %d hits from %s (%d passed quality filter)",
        len(hits), blast_tsv, len(hits),
    )
    return hits
