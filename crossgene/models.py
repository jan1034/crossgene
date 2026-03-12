"""Core data structures for crossgene."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class BedRegion:
    """A region from a BED file."""

    chrom: str
    start: int  # 0-based
    end: int  # exclusive
    name: str = "."  # column 4
    score: int = 0  # column 5
    strand: str = "."  # column 6


@dataclass
class GeneFeature:
    """A genomic feature (exon, CDS, UTR, etc.) within a gene."""

    feature_type: str  # "exon", "CDS", "UTR", etc.
    start: int  # genomic start (0-based)
    end: int  # genomic end


@dataclass
class GeneRecord:
    """A gene with its genomic coordinates, sequence, and features."""

    name: str  # gene name (e.g., "BRCA1")
    gene_id: str  # Ensembl ID (e.g., "ENSG00000012048")
    chrom: str  # chromosome (e.g., "chr17")
    start: int  # genomic start (0-based)
    end: int  # genomic end
    strand: str  # '+' or '-' (original strand from GTF)
    sequence: str = ""  # genomic-orientation sequence
    features: list[GeneFeature] = field(default_factory=list)
    gene_body_start: int = -1  # original gene start before flanking (-1 = not set)
    gene_body_end: int = -1  # original gene end before flanking (-1 = not set)


@dataclass
class AlignmentHit:
    """A single alignment hit between a query fragment and a target gene."""

    # Query gene coordinates (genomic)
    query_chrom: str
    query_start: int  # genomic start of fragment
    query_end: int  # genomic end of fragment

    # Target gene coordinates (genomic)
    target_chrom: str
    target_start: int  # hit start in target (genomic)
    target_end: int  # hit end in target (genomic)

    # Alignment info
    strand: str  # '+' or '-'
    identity: float  # 0.0 - 1.0
    query_coverage: float  # 0.0 - 1.0, fraction of query that aligned
    mapq: int  # mapping quality (pseudo-MAPQ from bitscore)
    is_primary: bool  # primary vs secondary alignment
    evalue: float = 0.0  # BLAST E-value
    bitscore: float = 0.0  # BLAST bitscore

    # Metadata
    query_gene: str = ""  # gene name of query
    target_gene: str = ""  # gene name of target
    direction: str = ""  # "A→B" or "B→A"


@dataclass
class FilterParams:
    """Filter thresholds for alignment hits."""

    min_quality: float = 30
    min_mapq: int = 0
    min_length: int = 25
    min_bitscore: float = 0.0
    max_evalue: float = float("inf")
