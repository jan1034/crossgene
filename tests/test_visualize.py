"""Tests for visualize module."""

import os

import pytest

from compare_genes.models import AlignmentHit, GeneFeature, GeneRecord
from compare_genes.visualize import (
    _mapq_to_alpha,
    _strand_label,
    _subsample_hits,
    create_circlize_plot,
)


def _make_gene(name, chrom, start, end, strand="+", features=None):
    return GeneRecord(
        name=name, gene_id="G001", chrom=chrom,
        start=start, end=end, strand=strand,
        sequence="A" * (end - start),
        features=features or [],
    )


def _make_hit(q_start, q_end, t_start, t_end, strand="+", mapq=40, score=90):
    return AlignmentHit(
        query_chrom="chr1", query_start=q_start, query_end=q_end,
        target_chrom="chr2", target_start=t_start, target_end=t_end,
        strand=strand, identity=0.9, mapq=mapq, cigar="50M",
        alignment_score=score, is_primary=True,
        query_gene="GENE_A", target_gene="GENE_B", direction="A→B",
    )


class TestHelpers:
    def test_strand_label_plus(self):
        gene = _make_gene("BRCA1", "chr17", 0, 100, "+")
        assert _strand_label(gene) == "BRCA1 (+)"

    def test_strand_label_minus(self):
        gene = _make_gene("TP53", "chr17", 0, 100, "-")
        assert _strand_label(gene) == "TP53 (-)"

    def test_mapq_to_alpha_range(self):
        assert _mapq_to_alpha(0, 60) == pytest.approx(0.2)
        assert _mapq_to_alpha(60, 60) == pytest.approx(1.0)
        assert _mapq_to_alpha(30, 60) == pytest.approx(0.6)

    def test_mapq_to_alpha_zero_max(self):
        assert _mapq_to_alpha(0, 0) == pytest.approx(0.6)

    def test_subsample(self):
        hits = [_make_hit(0, 50, 0, 50, score=i) for i in range(10)]
        result = _subsample_hits(hits, 3)
        assert len(result) == 3
        assert result[0].alignment_score == 9
        assert result[1].alignment_score == 8

    def test_subsample_under_limit(self):
        hits = [_make_hit(0, 50, 0, 50, score=i) for i in range(3)]
        result = _subsample_hits(hits, 5)
        assert len(result) == 3


class TestCreateCirclizePlot:
    def test_smoke_test(self, tmp_path):
        """Create a plot from synthetic data, verify PDF exists and is non-empty."""
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000, "+")
        gene_b = _make_gene("GENE_B", "chr2", 0, 800, "-")

        hits = [
            _make_hit(100, 150, 200, 250, strand="+", mapq=40, score=90),
            _make_hit(300, 350, 400, 450, strand="-", mapq=20, score=60),
            _make_hit(500, 550, 100, 150, strand="+", mapq=60, score=120),
        ]

        output = str(tmp_path / "test.pdf")
        create_circlize_plot(hits, gene_a, gene_b, output)

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_with_features(self, tmp_path):
        """Plot with gene features (exons)."""
        features_a = [
            GeneFeature("exon", 100, 200, {}),
            GeneFeature("CDS", 150, 200, {}),
        ]
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000, "+", features=features_a)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800, "-")

        hits = [_make_hit(100, 150, 200, 250)]
        output = str(tmp_path / "test.pdf")
        create_circlize_plot(hits, gene_a, gene_b, output)

        assert os.path.getsize(output) > 0

    def test_empty_hits(self, tmp_path):
        """Plot with no hits still produces a valid PDF."""
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800)

        output = str(tmp_path / "test.pdf")
        create_circlize_plot([], gene_a, gene_b, output)

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0
