"""PAF format parser with coordinate translation to genomic positions."""

from __future__ import annotations

import logging
from pathlib import Path

from crossgene.fragment import fragment_index_to_genomic
from crossgene.models import AlignmentHit, GeneRecord

logger = logging.getLogger(__name__)


def _parse_optional_fields(fields: list[str]) -> dict[str, str]:
    """Parse PAF optional fields (column 12+) into a dict.

    Optional fields have format TAG:TYPE:VALUE (e.g. cg:Z:50M, AS:i:100).
    """
    opts: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            opts[parts[0]] = parts[2]
    return opts


def _target_local_to_genomic(
    target_gene: GeneRecord, local_start: int, local_end: int
) -> tuple[int, int]:
    """Convert target-local coordinates to genomic coordinates.

    The target FASTA contains the sense sequence. For plus-strand genes,
    local position 0 corresponds to gene.start. For minus-strand genes,
    local position 0 corresponds to gene.end (the 5' end in sense
    orientation), so we reverse-map.
    """
    if target_gene.strand == "-":
        genomic_end = target_gene.end - local_start
        genomic_start = target_gene.end - local_end
    else:
        genomic_start = target_gene.start + local_start
        genomic_end = target_gene.start + local_end
    return genomic_start, genomic_end


def parse_paf(
    paf_path: Path,
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    fragment_size: int,
    step_size: int,
    min_quality: int,
    direction: str,
    min_mapq: int = 0,
) -> list[AlignmentHit]:
    """Parse a minimap2 PAF file into AlignmentHit objects.

    Translates fragment-local and target-local coordinates to genomic
    coordinates. Filters hits below min_quality (identity * 100).

    Args:
        paf_path: Path to PAF file.
        query_gene: The fragmented (query) gene with genomic coords.
        target_gene: The alignment reference gene with genomic coords.
        fragment_size: Fragment size used for fragmentation.
        step_size: Step size used for fragmentation.
        min_quality: Minimum identity threshold (0-100).
        direction: Label string, e.g. "A→B".

    Returns:
        List of AlignmentHit objects passing the quality filter.
    """
    hits: list[AlignmentHit] = []

    with open(paf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                logger.warning("Skipping malformed PAF line: %s", line[:80])
                continue

            # Core PAF columns
            query_name = fields[0]
            target_local_start = int(fields[7])
            target_local_end = int(fields[8])
            strand = fields[4]
            num_matches = int(fields[9])
            alignment_block_length = int(fields[10])
            mapq = int(fields[11])

            # Compute identity
            if alignment_block_length == 0:
                continue
            identity = num_matches / alignment_block_length

            # Filter by quality
            if identity * 100 < min_quality:
                continue

            # Filter by MAPQ
            if mapq < min_mapq:
                continue

            # Parse optional fields
            opts = _parse_optional_fields(fields[12:])
            cigar = opts.get("cg", "")
            alignment_score = int(opts["AS"]) if "AS" in opts else 0
            tp = opts.get("tp", "P")
            is_primary = tp == "P"

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
                target_gene, target_local_start, target_local_end
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
                    cigar=cigar,
                    alignment_score=alignment_score,
                    is_primary=is_primary,
                    query_gene=query_gene.name,
                    target_gene=target_gene.name,
                    direction=direction,
                )
            )

    logger.debug(
        "Parsed %d hits from %s (%d passed quality filter)",
        len(hits), paf_path, len(hits),
    )
    return hits
