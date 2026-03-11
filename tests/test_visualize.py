"""Tests for visualize module."""

import os

import pytest

from crossgene.models import AlignmentHit, GeneFeature, GeneRecord
from crossgene.models import BedRegion
from crossgene.visualize import (
    BedTrackConfig,
    _identity_to_alpha,
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
        strand=strand, identity=0.9, query_coverage=1.0, mapq=mapq, cigar="50M",
        alignment_score=score, is_primary=True,
        query_gene="GENE_A", target_gene="GENE_B", direction="A→B",
    )


class TestHelpers:
    @pytest.mark.parametrize("strand,expected", [("+", "BRCA1 (+)"), ("-", "BRCA1 (-)")])
    def test_strand_label(self, strand, expected):
        gene = _make_gene("BRCA1", "chr17", 0, 100, strand)
        assert _strand_label(gene) == expected

    def test_identity_to_alpha(self):
        # Few hits: density_factor=1.0
        assert _identity_to_alpha(0.0, 10) == pytest.approx(0.1)
        assert _identity_to_alpha(1.0, 10) == pytest.approx(1.0)
        assert _identity_to_alpha(0.5, 10) == pytest.approx(0.55)
        # Many hits: density_factor=0.1
        assert _identity_to_alpha(1.0, 1000) == pytest.approx(0.1)

    def test_subsample(self):
        hits = [_make_hit(0, 50, 0, 50, score=i) for i in range(10)]
        result = _subsample_hits(hits, 3)
        assert len(result) == 3
        assert result[0].alignment_score == 9  # top scores kept

        # Under limit: no change
        assert len(_subsample_hits(hits[:3], 5)) == 3


class TestCreateCirclizePlot:
    def test_smoke_test(self, tmp_path):
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000, "+")
        gene_b = _make_gene("GENE_B", "chr2", 0, 800, "-")
        hits = [
            _make_hit(100, 150, 200, 250, strand="+"),
            _make_hit(300, 350, 400, 450, strand="-"),
        ]
        output = str(tmp_path / "test.pdf")
        create_circlize_plot(hits, gene_a, gene_b, output)
        assert os.path.getsize(output) > 0

    def test_with_features(self, tmp_path):
        features_a = [GeneFeature("exon", 100, 200, {}), GeneFeature("CDS", 150, 200, {})]
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000, "+", features=features_a)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800, "-")
        output = str(tmp_path / "test.pdf")
        create_circlize_plot([_make_hit(100, 150, 200, 250)], gene_a, gene_b, output)
        assert os.path.getsize(output) > 0

    def test_empty_hits(self, tmp_path):
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800)
        output = str(tmp_path / "test.pdf")
        create_circlize_plot([], gene_a, gene_b, output)
        assert os.path.exists(output) and os.path.getsize(output) > 0


class TestLegend:
    def _get_legend_labels(self, tmp_path, gene_a, gene_b, hits, bed_tracks=None):
        import matplotlib.pyplot as plt
        output = str(tmp_path / "legend_test.pdf")
        create_circlize_plot(hits, gene_a, gene_b, output, bed_tracks=bed_tracks)
        fig = plt.gcf()
        return [t.get_text() for t in fig.legends[0].get_texts()]

    def test_legend_entries(self, tmp_path):
        """Legend has alignment entries; feature entries only when drawn."""
        features_a = [GeneFeature("exon", 100, 200, {})]
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000, "+", features=features_a)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800, "-")
        hits = [_make_hit(100, 150, 200, 250)]

        labels = self._get_legend_labels(tmp_path, gene_a, gene_b, hits)
        assert "Same-sense alignment" in labels
        assert "Antisense alignment" in labels
        assert "exon" in labels
        assert "CDS" not in labels  # not drawn, so not in legend


class TestBedTracks:
    def test_plot_with_bed_tracks(self, tmp_path):
        gene_a = _make_gene("GENE_A", "chr1", 0, 1000)
        gene_b = _make_gene("GENE_B", "chr2", 0, 800)
        hits = [_make_hit(100, 150, 200, 250)]
        bt = BedTrackConfig(
            label="my_repeats", color="darkorchid",
            regions_a=[BedRegion("chr1", 100, 300, "AluSc")],
            regions_b=[BedRegion("chr2", 200, 400, "L1")],
        )
        output = str(tmp_path / "test.pdf")
        create_circlize_plot(hits, gene_a, gene_b, output, bed_tracks=[bt])
        assert os.path.getsize(output) > 0

        # Verify legend includes bed track
        import matplotlib.pyplot as plt
        labels = [t.get_text() for t in plt.gcf().legends[0].get_texts()]
        assert "my_repeats" in labels
