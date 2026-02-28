"""Circlize visualization of gene comparison results."""

from __future__ import annotations

import logging

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for PDF output

from pycirclize import Circos

from compare_genes.models import AlignmentHit, GeneRecord

logger = logging.getLogger(__name__)

# Colors
COLOR_SAME_SENSE = "steelblue"
COLOR_ANTISENSE = "firebrick"
COLOR_EXON = "orange"
COLOR_CDS = "forestgreen"


def _strand_label(gene: GeneRecord) -> str:
    """Format gene label with strand: 'BRCA1 (+)' or 'TP53 (-)'."""
    return f"{gene.name} ({gene.strand})"


def _mapq_to_alpha(mapq: int, max_mapq: int) -> float:
    """Normalize mapq to alpha range [0.2, 1.0]."""
    if max_mapq <= 0:
        return 0.6
    return 0.2 + 0.8 * min(mapq, max_mapq) / max_mapq


def _subsample_hits(
    hits: list[AlignmentHit], max_arcs: int
) -> list[AlignmentHit]:
    """Keep top max_arcs hits by alignment_score if over limit."""
    if len(hits) <= max_arcs:
        return hits
    logger.info(
        "Subsampling arcs: %d -> %d (keeping top by alignment_score)",
        len(hits), max_arcs,
    )
    return sorted(hits, key=lambda h: h.alignment_score, reverse=True)[:max_arcs]


def create_circlize_plot(
    hits_ab: list[AlignmentHit],
    gene_a: GeneRecord,
    gene_b: GeneRecord,
    output_path: str,
    max_arcs: int = 5000,
) -> None:
    """Create a circular comparison plot and save as PDF.

    Gene A is displayed on the right, Gene B on the left.
    Arcs connect matching regions colored by strand (blue = same-sense,
    red = antisense) with alpha proportional to mapping quality.

    Args:
        hits_ab: Alignment hits from A→B direction.
        gene_a: Query gene record.
        gene_b: Target gene record.
        output_path: Output PDF path.
        max_arcs: Maximum arcs to draw (auto-subsample if exceeded).
    """
    hits_ab = _subsample_hits(hits_ab, max_arcs)

    gene_a_len = gene_a.end - gene_a.start
    gene_b_len = gene_b.end - gene_b.start

    label_a = _strand_label(gene_a)
    label_b = _strand_label(gene_b)

    # Create circos with two sectors; Gene A right (start=315), Gene B left
    circos = Circos(
        {label_a: gene_a_len, label_b: gene_b_len},
        space=15,
    )

    # Draw gene tracks with labels
    for sector in circos.sectors:
        # Outer track for gene body
        track = sector.add_track((95, 100))
        track.axis(fc="lightgrey", ec="black", lw=0.5)

        # Determine which gene this sector represents
        if sector.name == label_a:
            gene = gene_a
        else:
            gene = gene_b

        # Draw feature annotations if available
        for feat in gene.features:
            local_start = feat.start - gene.start
            local_end = feat.end - gene.start
            fc = COLOR_CDS if feat.feature_type == "CDS" else COLOR_EXON
            track.rect(local_start, local_end, fc=fc, ec="none", lw=0)

        # Sector label
        sector.text(sector.name, r=115, size=10, fontweight="bold")

        # Inner tick track
        inner_track = sector.add_track((90, 95))
        gene_len = gene.end - gene.start
        # Choose tick interval based on gene length
        if gene_len > 100_000:
            tick_interval = 20_000
        elif gene_len > 10_000:
            tick_interval = 5_000
        else:
            tick_interval = 1_000
        inner_track.xticks_by_interval(
            tick_interval,
            label_formatter=lambda v: f"{v/1000:.0f}kb",
            label_size=6,
        )

    # Find max mapq for alpha normalization
    max_mapq = max((h.mapq for h in hits_ab), default=1)

    # Draw arcs (links) between matching regions
    for hit in hits_ab:
        # Convert genomic coords to sector-local coords
        a_local_start = hit.query_start - gene_a.start
        a_local_end = hit.query_end - gene_a.start
        b_local_start = hit.target_start - gene_b.start
        b_local_end = hit.target_end - gene_b.start

        # Clamp to sector bounds
        a_local_start = max(0, min(a_local_start, gene_a_len))
        a_local_end = max(0, min(a_local_end, gene_a_len))
        b_local_start = max(0, min(b_local_start, gene_b_len))
        b_local_end = max(0, min(b_local_end, gene_b_len))

        color = COLOR_ANTISENSE if hit.strand == "-" else COLOR_SAME_SENSE
        alpha = _mapq_to_alpha(hit.mapq, max_mapq)

        circos.link(
            (label_a, a_local_start, a_local_end),
            (label_b, b_local_start, b_local_end),
            color=color,
            alpha=alpha,
        )

    fig = circos.plotfig()
    fig.savefig(output_path, bbox_inches="tight")
    logger.info("Saved circlize plot to %s", output_path)
