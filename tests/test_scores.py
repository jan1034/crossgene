"""Tests for scores module."""

import numpy as np
import pytest

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
        strand="+", mapq=40, cigar="50M",
        alignment_score=90, is_primary=True,
        query_gene="TEST", target_gene="OTHER", direction="A→B",
    )
    defaults.update(kwargs)
    return AlignmentHit(
        query_start=q_start, query_end=q_end,
        identity=identity, **defaults,
    )


class TestComputeScores:
    def test_single_fragment_perfect(self):
        """One fragment with identity=1.0 → covered bases get score 100."""
        gene = _make_gene(start=1000, end=1100)
        # Fragment covering genomic [1000, 1050)
        hits = [_make_hit(1000, 1050, 1.0)]

        scores = compute_scores(hits, gene, 50, 10)

        assert scores.shape == (100,)
        np.testing.assert_array_equal(scores[:50], 100.0)
        np.testing.assert_array_equal(scores[50:], 0.0)

    def test_overlapping_fragments_average(self):
        """Two overlapping fragments → bases in overlap get the mean."""
        gene = _make_gene(start=0, end=100)
        # Fragment 0: [0, 50), identity=1.0
        # Fragment 1: [10, 60), identity=0.5
        hits = [
            _make_hit(0, 50, 1.0),
            _make_hit(10, 60, 0.5),
        ]

        scores = compute_scores(hits, gene, 50, 10)

        # [0,10): covered by frag 0 only → 100
        np.testing.assert_allclose(scores[:10], 100.0)
        # [10,50): covered by both → mean(1.0, 0.5) * 100 = 75
        np.testing.assert_allclose(scores[10:50], 75.0)
        # [50,60): covered by frag 1 only → 50
        np.testing.assert_allclose(scores[50:60], 50.0)
        # [60,100): uncovered → 0
        np.testing.assert_array_equal(scores[60:], 0.0)

    def test_uncovered_bases_zero(self):
        """Bases not covered by any fragment get score 0."""
        gene = _make_gene(start=0, end=100)
        hits = []
        scores = compute_scores(hits, gene, 50, 10)
        np.testing.assert_array_equal(scores, 0.0)

    def test_best_identity_per_fragment(self):
        """Multiple hits for same fragment → only best identity is used."""
        gene = _make_gene(start=0, end=50)
        # Two hits for the same fragment [0, 50): identities 0.8 and 0.95
        hits = [
            _make_hit(0, 50, 0.8),
            _make_hit(0, 50, 0.95),
        ]

        scores = compute_scores(hits, gene, 50, 10)

        # Best identity = 0.95 → score = 95
        np.testing.assert_allclose(scores, 95.0)

    def test_non_overlapping_fragments(self):
        """Two non-overlapping fragments with different identities."""
        gene = _make_gene(start=0, end=100)
        hits = [
            _make_hit(0, 50, 0.9),
            _make_hit(50, 100, 0.6),
        ]

        scores = compute_scores(hits, gene, 50, 50)

        np.testing.assert_allclose(scores[:50], 90.0)
        np.testing.assert_allclose(scores[50:], 60.0)

    def test_three_overlapping_fragments(self):
        """Three fragments overlapping at one position."""
        gene = _make_gene(start=0, end=50)
        # All three cover [0, 50), identities 1.0, 0.8, 0.6
        # But they should have different query ranges to be different fragments
        hits = [
            _make_hit(0, 30, 1.0),   # frag covering [0,30)
            _make_hit(10, 40, 0.8),  # frag covering [10,40)
            _make_hit(20, 50, 0.6),  # frag covering [20,50)
        ]

        scores = compute_scores(hits, gene, 30, 10)

        # [0,10): only frag 0 → 100
        np.testing.assert_allclose(scores[:10], 100.0)
        # [10,20): frag 0 + frag 1 → mean(1.0, 0.8)*100 = 90
        np.testing.assert_allclose(scores[10:20], 90.0)
        # [20,30): all three → mean(1.0, 0.8, 0.6)*100 = 80
        np.testing.assert_allclose(scores[20:30], 80.0)
        # [30,40): frag 1 + frag 2 → mean(0.8, 0.6)*100 = 70
        np.testing.assert_allclose(scores[30:40], 70.0)
        # [40,50): frag 2 only → 60
        np.testing.assert_allclose(scores[40:50], 60.0)
