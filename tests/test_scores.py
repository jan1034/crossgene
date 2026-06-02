"""Tests for scores module (HSP-based scoring)."""

import numpy as np

from crossgene.models import AlignmentHit, GeneRecord
from crossgene.scores import compute_scores


def _make_gene(start=1000, end=1100, strand="+") -> GeneRecord:
    return GeneRecord(
        name="TEST", gene_id="G001", chrom="chr1",
        start=start, end=end, strand=strand,
        sequence="A" * (end - start),
    )


def _make_hit(q_start, q_end, identity, **kwargs) -> AlignmentHit:
    defaults = dict(
        query_chrom="chr1", target_chrom="chr2",
        target_start=5000, target_end=5050,
        strand="+", query_coverage=1.0, mapq=40,
        is_primary=True,
        query_gene="TEST", target_gene="OTHER", direction="A→B",
    )
    defaults.update(kwargs)
    return AlignmentHit(
        query_start=q_start, query_end=q_end,
        identity=identity, **defaults,
    )


class TestComputeScores:
    def test_single_hsp_perfect(self):
        gene = _make_gene(start=1000, end=1100)
        scores = compute_scores([_make_hit(1000, 1050, 1.0)], gene)
        assert scores.shape == (100,)
        np.testing.assert_array_equal(scores[:50], 100.0)
        np.testing.assert_array_equal(scores[50:], 0.0)

    def test_overlapping_hsps_best_wins(self):
        """Overlapping HSPs: max identity per base, not average."""
        gene = _make_gene(start=0, end=100)
        hits = [_make_hit(0, 50, 1.0), _make_hit(0, 50, 0.8), _make_hit(10, 60, 0.5)]
        scores = compute_scores(hits, gene)

        np.testing.assert_allclose(scores[:10], 100.0)   # max(1.0, 0.8) = 1.0
        np.testing.assert_allclose(scores[10:50], 100.0)  # max(1.0, 0.8, 0.5) = 1.0
        np.testing.assert_allclose(scores[50:60], 50.0)   # only 0.5 HSP
        np.testing.assert_array_equal(scores[60:], 0.0)

    def test_uncovered_bases_zero(self):
        gene = _make_gene(start=0, end=100)
        scores = compute_scores([], gene)
        np.testing.assert_array_equal(scores, 0.0)

    def test_non_overlapping_hsps(self):
        gene = _make_gene(start=0, end=100)
        hits = [_make_hit(0, 50, 0.9), _make_hit(50, 100, 0.6)]
        scores = compute_scores(hits, gene)
        np.testing.assert_allclose(scores[:50], 90.0)
        np.testing.assert_allclose(scores[50:], 60.0)

    def test_three_overlapping_hsps_max(self):
        """Three HSPs: max identity per base position."""
        gene = _make_gene(start=0, end=50)
        hits = [_make_hit(0, 30, 1.0), _make_hit(10, 40, 0.8), _make_hit(20, 50, 0.6)]
        scores = compute_scores(hits, gene)

        np.testing.assert_allclose(scores[:20], 100.0)    # max includes 1.0
        np.testing.assert_allclose(scores[20:30], 100.0)  # all three, max=1.0
        np.testing.assert_allclose(scores[30:40], 80.0)   # max(0.8, 0.6)
        np.testing.assert_allclose(scores[40:50], 60.0)   # only 0.6
