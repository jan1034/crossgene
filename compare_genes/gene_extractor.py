"""Gene lookup in GTF and sequence extraction from genome FASTA."""

from __future__ import annotations

import logging

import gtfparse
import pysam

from compare_genes.models import GeneFeature, GeneRecord

logger = logging.getLogger(__name__)


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def lookup_gene(gene_name: str, gtf_path: str) -> GeneRecord:
    """Look up a gene by name in a GTF file.

    Searches the gene_name attribute. If multiple matches are found,
    logs a warning and returns the first. Raises ValueError if not found.

    Returns a GeneRecord without sequence (sequence must be filled via
    extract_sequence).

    Note: GTF coordinates are 1-based inclusive. We convert to 0-based
    half-open internally (start = gtf_start - 1, end = gtf_end).
    """
    df = gtfparse.read_gtf(gtf_path, result_type="pandas")

    # Filter for matching gene_name
    genes = df[(df["feature"] == "gene") & (df["gene_name"] == gene_name)]

    if len(genes) == 0:
        raise ValueError(
            f"Gene '{gene_name}' not found in GTF. "
            f"Check spelling (search is case-sensitive)."
        )

    if len(genes) > 1:
        locations = [
            f"  {row['seqname']}:{row['start']}-{row['end']} ({row['gene_id']})"
            for _, row in genes.iterrows()
        ]
        logger.warning(
            "Multiple matches for gene '%s', using first:\n%s",
            gene_name,
            "\n".join(locations),
        )

    gene_row = genes.iloc[0]

    # Convert GTF 1-based inclusive to 0-based half-open
    start_0based = int(gene_row["start"]) - 1
    end_0based = int(gene_row["end"])

    # Extract sub-gene features (exon, CDS, UTR, etc.) if available
    features = []
    gene_id = str(gene_row["gene_id"])
    sub_features = df[
        (df["gene_id"] == gene_id) & (df["feature"] != "gene")
    ]
    for _, feat_row in sub_features.iterrows():
        feature_type = str(feat_row["feature"])
        if feature_type in ("exon", "CDS", "five_prime_utr", "three_prime_utr"):
            features.append(
                GeneFeature(
                    feature_type=feature_type,
                    start=int(feat_row["start"]) - 1,  # 0-based
                    end=int(feat_row["end"]),
                    metadata={},
                )
            )

    return GeneRecord(
        name=gene_name,
        gene_id=gene_id,
        chrom=str(gene_row["seqname"]),
        start=start_0based,
        end=end_0based,
        strand=str(gene_row["strand"]),
        sequence="",
        features=sorted(features, key=lambda f: f.start),
    )


def extract_sequence(gene: GeneRecord, genome_path: str) -> GeneRecord:
    """Extract the gene sequence from a genome FASTA.

    Fetches chrom:start-end using pysam (0-based half-open coords).
    For minus-strand genes, reverse-complements the sequence to produce
    the sense strand.

    Returns a new GeneRecord with the sequence populated.
    """
    fa = pysam.FastaFile(genome_path)
    try:
        seq = fa.fetch(gene.chrom, gene.start, gene.end)
    finally:
        fa.close()

    if gene.strand == "-":
        seq = _reverse_complement(seq)

    return GeneRecord(
        name=gene.name,
        gene_id=gene.gene_id,
        chrom=gene.chrom,
        start=gene.start,
        end=gene.end,
        strand=gene.strand,
        sequence=seq,
        features=gene.features,
    )
