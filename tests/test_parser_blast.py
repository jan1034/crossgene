"""Tests for parser_blast module (gene-vs-gene)."""

import pytest

from crossgene.models import GeneRecord
from crossgene.parser_blast import parse_blast_gene_vs_gene


def _make_gene(name="GENE_A", chrom="chr1", start=1000, end=2000, strand="+") -> GeneRecord:
    return GeneRecord(
        name=name, gene_id="G001", chrom=chrom,
        start=start, end=end, strand=strand,
        sequence="A" * (end - start),
    )


def _blast_line(
    qseqid="GENE_A", sseqid="GENE_B", pident="90.000", length="200",
    mismatch="20", gapopen="0", qstart="51", qend="250",
    sstart="101", send="300", evalue="1e-40", bitscore="180.0",
    qlen="1000", slen="1000", sstrand="plus",
):
    return "\t".join([
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qlen, slen, sstrand,
    ])


class TestParseBasic:
    def test_single_hsp(self, tmp_path):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line() + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, 0, "A→B")

        assert len(hits) == 1
        h = hits[0]
        assert h.query_gene == "GENE_A"
        assert h.target_gene == "GENE_B"
        assert h.identity == 0.9
        assert h.is_primary is True
        assert h.evalue == pytest.approx(1e-40)
        assert h.bitscore == 180.0

    def test_multiple_hsps_primary_secondary(self, tmp_path):
        """Multiple HSPs: first is primary, rest secondary."""
        tsv = tmp_path / "test.tsv"
        lines = [
            _blast_line(qstart="1", qend="200", bitscore="90.0"),
            _blast_line(qstart="300", qend="500", bitscore="70.0"),
            _blast_line(qstart="600", qend="800", bitscore="80.0"),
        ]
        tsv.write_text("\n".join(lines) + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, 0, "A→B")
        assert len(hits) == 3
        assert hits[0].is_primary is True
        assert hits[1].is_primary is False
        assert hits[2].is_primary is False


class TestFilters:
    @pytest.mark.parametrize("min_qual,filter_kwargs,expected_count", [
        (30, {}, 1),                      # filters pident=20
        (0, {"min_length": 25}, 1),       # filters length=20
        (0, {"min_bitscore": 50.0}, 1),   # filters bitscore=10
        (0, {"max_evalue": 0.01}, 1),     # filters evalue=0.5
    ])
    def test_filter(self, tmp_path, min_qual, filter_kwargs, expected_count):
        tsv = tmp_path / "test.tsv"
        good = _blast_line()
        bad = _blast_line(pident="20.000", length="20", bitscore="10.0", evalue="0.5", qstart="300", qend="319")
        tsv.write_text(good + "\n" + bad + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, min_qual, "A→B", **filter_kwargs)
        assert len(hits) == expected_count


class TestCoordinateTranslation:
    @pytest.mark.parametrize("gene_strand,qstart,qend,expected_start,expected_end", [
        ("+", "51", "250", 1050, 1250),   # plus strand query
        ("-", "51", "250", 1050, 1250),   # minus strand query (genomic orientation, same mapping)
    ])
    def test_query_coords(self, tmp_path, gene_strand, qstart, qend, expected_start, expected_end):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(qstart=qstart, qend=qend) + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, gene_strand)
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, 0, "A→B")
        assert hits[0].query_start == expected_start
        assert hits[0].query_end == expected_end

    @pytest.mark.parametrize("gene_strand,sstart,send,sstrand,expected_start,expected_end,expected_strand", [
        ("+", "101", "300", "plus", 5100, 5300, "+"),    # plus target, plus hit → same-sense
        ("-", "101", "300", "plus", 5100, 5300, "-"),    # minus target, plus hit → antisense (query+ target-)
        ("+", "300", "101", "minus", 5100, 5300, "-"),   # minus strand BLAST hit
    ])
    def test_target_coords(self, tmp_path, gene_strand, sstart, send, sstrand,
                           expected_start, expected_end, expected_strand):
        tsv = tmp_path / "test.tsv"
        tsv.write_text(_blast_line(sstart=sstart, send=send, sstrand=sstrand) + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, gene_strand)

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, 0, "A→B")
        assert hits[0].target_start == expected_start
        assert hits[0].target_end == expected_end
        assert hits[0].strand == expected_strand


class TestMapqFromBitscore:
    def test_mapq_scaling(self, tmp_path):
        """Max bitscore → MAPQ 60, half bitscore → MAPQ 30."""
        tsv = tmp_path / "test.tsv"
        lines = [
            _blast_line(bitscore="100.0", qstart="1", qend="200"),
            _blast_line(bitscore="50.0", qstart="300", qend="500"),
        ]
        tsv.write_text("\n".join(lines) + "\n")
        query_gene = _make_gene("GENE_A", "chr1", 1000, 2000, "+")
        target_gene = _make_gene("GENE_B", "chr2", 5000, 6000, "+")

        hits = parse_blast_gene_vs_gene(tsv, query_gene, target_gene, 0, "A→B")
        assert hits[0].mapq == 60
        assert hits[1].mapq == 30


class TestEmptyOutput:
    def test_empty_or_comments(self, tmp_path):
        for content in ["", "# comment line\n# another\n"]:
            tsv = tmp_path / "test.tsv"
            tsv.write_text(content)
            hits = parse_blast_gene_vs_gene(tsv, _make_gene(), _make_gene("GENE_B", "chr2", 5000, 6000), 0, "A→B")
            assert hits == []
