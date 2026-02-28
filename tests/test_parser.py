"""Tests for parser module."""

from pathlib import Path

import pytest

from compare_genes.models import AlignmentHit, GeneRecord
from compare_genes.parser import _target_local_to_genomic, parse_paf


def _make_gene(name="GENE_A", chrom="chr1", start=1000, end=2000, strand="+") -> GeneRecord:
    seq_len = end - start
    return GeneRecord(
        name=name, gene_id="G001", chrom=chrom,
        start=start, end=end, strand=strand,
        sequence="A" * seq_len,
    )


# PAF line template: 12 mandatory + optional fields
# query_name query_len q_start q_end strand target_name target_len t_start t_end matches block_len mapq [optional]
def _paf_line(
    query_name="GENE_A:0", query_len=50, q_start=0, q_end=50,
    strand="+", target_name="GENE_B", target_len=1000,
    t_start=100, t_end=150, matches=45, block_len=50, mapq=40,
    extra="cg:Z:50M\tAS:i:90\ttp:A:P",
):
    fields = [
        query_name, str(query_len), str(q_start), str(q_end),
        strand, target_name, str(target_len),
        str(t_start), str(t_end),
        str(matches), str(block_len), str(mapq),
    ]
    if extra:
        fields.append(extra)
    return "\t".join(fields)


class TestTargetLocalToGenomic:
    def test_plus_strand(self):
        gene = _make_gene(start=5000, end=6000, strand="+")
        assert _target_local_to_genomic(gene, 100, 150) == (5100, 5150)

    def test_minus_strand(self):
        gene = _make_gene(start=5000, end=6000, strand="-")
        # local [100, 150) → genomic [6000-150, 6000-100) = [5850, 5900)
        assert _target_local_to_genomic(gene, 100, 150) == (5850, 5900)


class TestParsePaf:
    def test_basic_parse(self, tmp_path):
        """Parse a single well-formed PAF line."""
        paf = tmp_path / "test.paf"
        paf.write_text(_paf_line() + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")

        assert len(hits) == 1
        h = hits[0]
        assert h.query_gene == "GENE_A"
        assert h.target_gene == "GENE_B"
        assert h.direction == "A→B"
        assert h.strand == "+"
        assert h.identity == 45 / 50  # 0.9
        assert h.mapq == 40
        assert h.cigar == "50M"
        assert h.alignment_score == 90
        assert h.is_primary is True

    def test_query_genomic_coords_plus_strand(self, tmp_path):
        """Fragment index 2, step 10 → sense [20,70) → genomic [1020, 1070)."""
        paf = tmp_path / "test.paf"
        paf.write_text(_paf_line(query_name="GENE_A:2") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].query_start == 1020
        assert hits[0].query_end == 1070

    def test_query_genomic_coords_minus_strand(self, tmp_path):
        """Minus-strand query: fragment 0 → genomic [end-50, end)."""
        paf = tmp_path / "test.paf"
        paf.write_text(_paf_line(query_name="GENE_A:0") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "-")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        # frag 0: sense [0,50) → minus strand genomic [2000-50, 2000) = [1950, 2000)
        assert hits[0].query_start == 1950
        assert hits[0].query_end == 2000

    def test_target_genomic_coords(self, tmp_path):
        """Target local [100, 150) with plus-strand target starting at 5000."""
        paf = tmp_path / "test.paf"
        paf.write_text(_paf_line(t_start=100, t_end=150) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].target_start == 5100
        assert hits[0].target_end == 5150

    def test_target_minus_strand(self, tmp_path):
        """Minus-strand target: local [100,150) → genomic [5850, 5900)."""
        paf = tmp_path / "test.paf"
        paf.write_text(_paf_line(t_start=100, t_end=150) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "-")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].target_start == 5850
        assert hits[0].target_end == 5900

    def test_quality_filter(self, tmp_path):
        """Hits below min_quality are excluded."""
        paf = tmp_path / "test.paf"
        # identity = 45/50 = 90%
        high_q = _paf_line(query_name="GENE_A:0", matches=45, block_len=50)
        # identity = 10/50 = 20%
        low_q = _paf_line(query_name="GENE_A:1", matches=10, block_len=50)
        paf.write_text(high_q + "\n" + low_q + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 30, "A→B")
        assert len(hits) == 1
        assert hits[0].identity == 0.9

    def test_secondary_alignment(self, tmp_path):
        """Secondary alignments have is_primary=False."""
        paf = tmp_path / "test.paf"
        line = _paf_line(extra="cg:Z:50M\tAS:i:90\ttp:A:S")
        paf.write_text(line + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].is_primary is False

    def test_multiple_lines(self, tmp_path):
        """Parse multiple PAF lines."""
        paf = tmp_path / "test.paf"
        lines = [
            _paf_line(query_name="GENE_A:0"),
            _paf_line(query_name="GENE_A:1", t_start=200, t_end=250),
            _paf_line(query_name="GENE_A:2", t_start=300, t_end=350),
        ]
        paf.write_text("\n".join(lines) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert len(hits) == 3

    def test_empty_paf(self, tmp_path):
        paf = tmp_path / "test.paf"
        paf.write_text("")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert len(hits) == 0

    def test_missing_optional_fields(self, tmp_path):
        """PAF line with no optional fields → defaults for cigar, score, primary."""
        paf = tmp_path / "test.paf"
        line = _paf_line(extra="")
        paf.write_text(line + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_paf(paf, query_gene, target_gene, 50, 10, 0, "A→B")
        assert len(hits) == 1
        assert hits[0].cigar == ""
        assert hits[0].alignment_score == 0
        assert hits[0].is_primary is True  # default when tp tag missing
