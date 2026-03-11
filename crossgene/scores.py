"""Per-base similarity scoring from BLAST HSPs."""

from __future__ import annotations

import numpy as np

from crossgene.models import AlignmentHit, GeneRecord


def compute_scores(
    hits: list[AlignmentHit],
    gene: GeneRecord,
) -> np.ndarray:
    """Compute per-base similarity score from BLAST HSPs.

    For each base in the query gene, assign the identity of the best
    overlapping HSP (max, not average). Scale to 0-100.
    Uncovered bases get score 0.

    Returns a numpy array of length (gene.end - gene.start) with scores 0-100.
    """
    gene_len = gene.end - gene.start
    scores = np.zeros(gene_len, dtype=np.float64)

    for hit in hits:
        start = hit.query_start - gene.start
        end = hit.query_end - gene.start

        # Clamp to array bounds
        start = max(0, start)
        end = min(gene_len, end)

        if start >= end:
            continue

        identity_score = hit.identity * 100.0
        mask = scores[start:end] < identity_score
        scores[start:end] = np.where(mask, identity_score, scores[start:end])

    return scores
