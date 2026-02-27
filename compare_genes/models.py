"""Core data structures for compare_genes."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class GeneFeature:
    """A genomic feature (exon, CDS, UTR, etc.) within a gene."""

    feature_type: str  # "exon", "CDS", "UTR", etc.
    start: int  # genomic start (0-based)
    end: int  # genomic end
    metadata: dict = field(default_factory=dict)


@dataclass
class GeneRecord:
    """A gene with its genomic coordinates, sequence, and features."""

    name: str  # gene name (e.g., "BRCA1")
    gene_id: str  # Ensembl ID (e.g., "ENSG00000012048")
    chrom: str  # chromosome (e.g., "chr17")
    start: int  # genomic start (0-based)
    end: int  # genomic end
    strand: str  # '+' or '-' (original strand from GTF)
    sequence: str = ""  # sense sequence (revcomp applied if minus strand)
    features: list[GeneFeature] = field(default_factory=list)


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
    mapq: int  # mapping quality
    cigar: str  # CIGAR string
    alignment_score: int  # AS tag value
    is_primary: bool  # primary vs secondary alignment

    # Metadata
    query_gene: str  # gene name of query
    target_gene: str  # gene name of target
    direction: str  # "A→B" or "B→A"
