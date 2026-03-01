"""Tests for parser_blast module."""

from pathlib import Path

import pytest

from crossgene.models import AlignmentHit, GeneRecord
from crossgene.parser_blast import parse_blast_tabular


def _make_gene(name="GENE_A", chrom="chr1", start=1000, end=2000, strand="+") -> GeneRecord:
    seq_len = end - start
    return GeneRecord(
        name=name, gene_id="G001", chrom=chrom,
        start=start, end=end, strand=strand,
        sequence="A" * seq_len,
    )


# BLAST tabular line: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand
def _blast_line(
    qseqid="GENE_A:0", sseqid="GENE_B", pident="90.000", length="50",
    mismatch="5", gapopen="0", qstart="1", qend="50",
    sstart="101", send="150", evalue="1e-10", bitscore="80.0",
    qlen="50", slen="1000", sstrand="plus",
):
    return "\t".join([
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qlen, slen, sstrand,
    ])


class TestParseBasic:
    def test_single_hit(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line() + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")

        assert len(hits) == 1
        h = hits[0]
        assert h.query_gene == "GENE_A"
        assert h.target_gene == "GENE_B"
        assert h.direction == "A→B"
        assert h.identity == 0.9
        assert h.cigar == ""
        assert h.is_primary is True

    def test_multiple_hits(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        lines = [
            _blast_line(qseqid="GENE_A:0"),
            _blast_line(qseqid="GENE_A:1", sstart="201", send="250"),
            _blast_line(qseqid="GENE_A:2", sstart="301", send="350"),
        ]
        tsv.write_text("\n".join(lines) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert len(hits) == 3


class TestIdentityFilter:
    def test_filters_low_identity(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        high_q = _blast_line(qseqid="GENE_A:0", pident="90.000")
        low_q = _blast_line(qseqid="GENE_A:1", pident="20.000")
        tsv.write_text(high_q + "\n" + low_q + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 30, "A→B")
        assert len(hits) == 1
        assert hits[0].identity == 0.9


class TestCoordinateTranslation:
    def test_plus_strand_target(self, tmp_path):
        """sstart=101, send=150 (1-based) → local [100, 150) → genomic [5100, 5150)."""
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(sstart="101", send="150") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].target_start == 5100
        assert hits[0].target_end == 5150

    def test_minus_strand_target(self, tmp_path):
        """Minus-strand target: local [100, 150) → genomic [5850, 5900)."""
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(sstart="101", send="150") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "-")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].target_start == 5850
        assert hits[0].target_end == 5900

    def test_query_coords_plus_strand(self, tmp_path):
        """Fragment index 2, step 10 → sense [20,70) → genomic [1020, 1070)."""
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(qseqid="GENE_A:2") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].query_start == 1020
        assert hits[0].query_end == 1070

    def test_query_coords_minus_strand(self, tmp_path):
        """Minus-strand query: fragment 0 → genomic [1950, 2000)."""
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(qseqid="GENE_A:0") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "-")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].query_start == 1950
        assert hits[0].query_end == 2000

    def test_minus_strand_blast_coords(self, tmp_path):
        """BLAST minus strand: sstart > send → swap and convert."""
        tsv = tmp_path / "test.tsv"
        # sstart=150 > send=101 indicates minus strand hit
        tsv.write_text(_blast_line(sstart="150", send="101", sstrand="minus") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].target_start == 5100
        assert hits[0].target_end == 5150
        assert hits[0].strand == "-"


class TestMapqFromBitscore:
    def test_max_bitscore_gets_60(self, tmp_path):
        """Single hit: bitscore is max, so mapq = 60."""
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(bitscore="100.0") + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].mapq == 60

    def test_half_bitscore(self, tmp_path):
        """Two hits: one with max bitscore (60), other with half (30)."""
        tsv = tmp_path / "test.tsv"
        lines = [
            _blast_line(qseqid="GENE_A:0", bitscore="100.0"),
            _blast_line(qseqid="GENE_A:1", bitscore="50.0"),
        ]
        tsv.write_text("\n".join(lines) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].mapq == 60
        assert hits[1].mapq == 30


class TestStrandMapping:
    def test_plus_strand(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(sstrand="plus") + "\n")

        query_gene = _make_gene()
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].strand == "+"

    def test_minus_strand(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(sstrand="minus", sstart="150", send="101") + "\n")

        query_gene = _make_gene()
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits[0].strand == "-"


class TestPrimarySecondary:
    def test_first_is_primary_rest_secondary(self, tmp_path):
        """First HSP per query is primary, subsequent are secondary."""
        tsv = tmp_path / "test.tsv"
        lines = [
            _blast_line(qseqid="GENE_A:0", sstart="101", send="150", bitscore="90.0"),
            _blast_line(qseqid="GENE_A:0", sstart="201", send="250", bitscore="70.0"),
            _blast_line(qseqid="GENE_A:1", sstart="301", send="350", bitscore="80.0"),
        ]
        tsv.write_text("\n".join(lines) + "\n")

        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert len(hits) == 3
        assert hits[0].is_primary is True   # first hit for GENE_A:0
        assert hits[1].is_primary is False  # second hit for GENE_A:0
        assert hits[2].is_primary is True   # first hit for GENE_A:1


class TestEmptyOutput:
    def test_empty_file(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text("")

        query_gene = _make_gene()
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits == []

    def test_only_comments(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text("# comment line\n# another comment\n")

        query_gene = _make_gene()
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_tabular(tsv, query_gene, target_gene, 50, 10, 0, "A→B")
        assert hits == []
