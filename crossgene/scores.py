"""Per-base mappability score aggregation."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from crossgene.models import AlignmentHit, GeneRecord


def compute_scores(
    hits: list[AlignmentHit],
    gene: GeneRecord,
    fragment_size: int,
    step_size: int,
) -> np.ndarray:
    """Compute per-base mappability scores (0-100) for a gene.

    Algorithm:
    1. Group hits by fragment (identified by query_start, query_end).
    2. For each fragment, take the best (highest) identity.
    3. For each genomic base covered by the fragment, accumulate that
       best identity and a coverage count.
    4. Final score per base = mean of best identities * 100.
    5. Uncovered bases get score 0.

    Returns a numpy array of length (gene.end - gene.start) with scores 0-100.
    """
    gene_len = gene.end - gene.start
    identity_sum = np.zeros(gene_len, dtype=np.float64)
    coverage = np.zeros(gene_len, dtype=np.float64)

    # Group hits by fragment (same query genomic range = same fragment)
    best_per_fragment: dict[tuple[int, int], float] = defaultdict(float)
    for hit in hits:
        key = (hit.query_start, hit.query_end)
        if hit.identity > best_per_fragment[key]:
            best_per_fragment[key] = hit.identity

    # Accumulate scores
    for (q_start, q_end), best_identity in best_per_fragment.items():
        arr_start = q_start - gene.start
        arr_end = q_end - gene.start

        # Clamp to array bounds (shouldn't be needed, but safe)
        arr_start = max(0, arr_start)
        arr_end = min(gene_len, arr_end)

        identity_sum[arr_start:arr_end] += best_identity
        coverage[arr_start:arr_end] += 1.0

    # Compute mean score, avoiding division by zero
    scores = np.zeros(gene_len, dtype=np.float64)
    mask = coverage > 0
    scores[mask] = (identity_sum[mask] / coverage[mask]) * 100.0

    return scores
