"""BLASTN tabular output parser for gene-vs-gene alignment."""

from __future__ import annotations

import logging
from pathlib import Path

from crossgene.models import AlignmentHit, GeneRecord

logger = logging.getLogger(__name__)

# Column indices for -outfmt "6 qseqid sseqid pident length mismatch gapopen
#                               qstart qend sstart send evalue bitscore qlen slen sstrand"
_COL_QSEQID = 0
_COL_PIDENT = 2
_COL_LENGTH = 3
_COL_QSTART = 6
_COL_QEND = 7
_COL_SSTART = 8
_COL_SEND = 9
_COL_EVALUE = 10
_COL_BITSCORE = 11
_COL_QLEN = 12
_COL_SSTRAND = 14
_NUM_COLUMNS = 15


def _local_to_genomic(
    gene: GeneRecord, local_start: int, local_end: int
) -> tuple[int, int]:
    """Convert sequence-local coordinates to genomic coordinates.

    The FASTA contains the genomic-orientation sequence, so local
    position 0 always corresponds to gene.start regardless of strand.
    """
    genomic_start = gene.start + local_start
    genomic_end = gene.start + local_end
    return genomic_start, genomic_end


def parse_blast_gene_vs_gene(
    blast_tsv: Path,
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    min_quality: float,
    direction: str,
    min_mapq: int = 0,
    min_length: int = 0,
    min_bitscore: float = 0.0,
    max_evalue: float = float("inf"),
) -> list[AlignmentHit]:
    """Parse BLASTN tabular output from gene-vs-gene alignment.

    Coordinate translation:
    - Query: BLAST qstart/qend (1-based inclusive) → 0-based half-open → genomic
    - Target: BLAST sstart/send (1-based inclusive) → 0-based half-open → genomic

    Args:
        blast_tsv: Path to BLASTN tabular output file.
        query_gene: The query gene with genomic coords.
        target_gene: The subject gene with genomic coords.
        min_quality: Minimum identity threshold (0-100).
        direction: Label string, e.g. "A→B".
        min_mapq: Minimum pseudo-MAPQ threshold.
        min_length: Minimum HSP alignment length.
        min_bitscore: Minimum bitscore threshold.
        max_evalue: Maximum E-value threshold.

    Returns:
        List of AlignmentHit objects passing all filters.
    """
    hits: list[AlignmentHit] = []
    max_bs = 0.0
    hsp_count = 0

    with open(blast_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < _NUM_COLUMNS:
                logger.warning("Skipping malformed BLAST line: %s", line[:80])
                continue

            pident = float(fields[_COL_PIDENT])  # 0-100 from BLAST
            identity = pident / 100.0
            alignment_length = int(fields[_COL_LENGTH])
            evalue = float(fields[_COL_EVALUE])
            bitscore = float(fields[_COL_BITSCORE])
            qlen = int(fields[_COL_QLEN])
            query_coverage = alignment_length / qlen if qlen > 0 else 0.0

            # Apply filters (except min_mapq, deferred to post-normalization)
            if pident < min_quality:
                continue
            if alignment_length < min_length:
                continue
            if bitscore < min_bitscore:
                continue
            if evalue > max_evalue:
                continue

            if bitscore > max_bs:
                max_bs = bitscore

            # Query coords: BLAST qstart/qend are 1-based inclusive, always qstart < qend
            qstart = int(fields[_COL_QSTART])
            qend = int(fields[_COL_QEND])
            q_local_start = qstart - 1  # 1-based → 0-based
            q_local_end = qend           # inclusive → exclusive

            # Target coords: BLAST sstart/send are 1-based inclusive
            # sstart > send for minus strand hits
            sstart = int(fields[_COL_SSTART])
            send = int(fields[_COL_SEND])
            sstrand = fields[_COL_SSTRAND]
            blast_strand = "+" if sstrand == "plus" else "-"

            # Derive relative strand (same-sense vs antisense).
            # Both sequences are genomic orientation. BLAST "plus" means same
            # genomic direction. Same-sense requires considering both genes' strands.
            if query_gene.strand == target_gene.strand:
                strand = blast_strand
            else:
                strand = "-" if blast_strand == "+" else "+"

            if sstart <= send:
                s_local_start = sstart - 1
                s_local_end = send
            else:
                s_local_start = send - 1
                s_local_end = sstart

            # Primary: first HSP is primary, rest are secondary
            is_primary = hsp_count == 0
            hsp_count += 1

            # Translate to genomic coordinates
            query_genomic_start, query_genomic_end = _local_to_genomic(
                query_gene, q_local_start, q_local_end
            )
            target_genomic_start, target_genomic_end = _local_to_genomic(
                target_gene, s_local_start, s_local_end
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
                    query_coverage=query_coverage,
                    mapq=0,  # placeholder, normalized below
                    is_primary=is_primary,
                    evalue=evalue,
                    bitscore=bitscore,
                    query_gene=query_gene.name,
                    target_gene=target_gene.name,
                    direction=direction,
                )
            )

    # Normalize MAPQ and apply min_mapq filter
    if max_bs == 0:
        max_bs = 1.0
    for h in hits:
        h.mapq = min(60, int(h.bitscore / max_bs * 60))
    if min_mapq > 0:
        hits = [h for h in hits if h.mapq >= min_mapq]

    if not hits:
        logger.debug("Empty BLAST output: %s", blast_tsv)

    logger.debug(
        "Parsed %d HSPs from %s (passing all filters)",
        len(hits), blast_tsv,
    )
    return hits
