"""Tests for tsv_writer module."""

import csv

from crossgene.models import AlignmentHit
from crossgene.tsv_writer import COLUMNS, write_tsv


def _make_hit(q_start, q_end, alignment_score, identity=0.9, **kwargs) -> AlignmentHit:
    defaults = dict(
        query_chrom="chr1", target_chrom="chr2",
        target_start=5000, target_end=5050,
        strand="+", query_coverage=1.0, mapq=40, cigar="50M",
        is_primary=True,
        query_gene="GENE_A", target_gene="GENE_B", direction="A→B",
    )
    defaults.update(kwargs)
    return AlignmentHit(
        query_start=q_start, query_end=q_end,
        identity=identity, alignment_score=alignment_score,
        **defaults,
    )


class TestWriteTsv:
    def test_header_and_columns(self, tmp_path):
        output = str(tmp_path / "test.tsv")
        write_tsv([], output)

        with open(output) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header == COLUMNS
        assert len(header) == 16

    def test_single_hit(self, tmp_path):
        output = str(tmp_path / "test.tsv")
        hits = [_make_hit(1000, 1050, 90, identity=0.92)]
        write_tsv(hits, output)

        with open(output) as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)  # skip header
            row = next(reader)

        assert row[0] == "GENE_A"      # query_gene
        assert row[1] == "chr1"        # query_chrom
        assert row[2] == "1000"        # frag_start
        assert row[3] == "1050"        # frag_end
        assert row[4] == "GENE_B"      # target_gene
        assert row[9] == "0.9200"      # identity
        assert row[10] == "1.0000"     # query_coverage
        assert row[13] == "true"       # is_primary
        assert row[15] == "A→B"        # direction

    def test_sorting(self, tmp_path):
        """Sorted by frag_start asc, then alignment_score desc."""
        output = str(tmp_path / "test.tsv")
        hits = [
            _make_hit(2000, 2050, 50),   # second position, low score
            _make_hit(1000, 1050, 80),   # first position, lower score
            _make_hit(1000, 1050, 100),  # first position, higher score
            _make_hit(2000, 2050, 90),   # second position, high score
        ]
        write_tsv(hits, output)

        with open(output) as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)
            rows = list(reader)

        assert len(rows) == 4
        # First position, highest score first
        assert rows[0][2] == "1000" and rows[0][12] == "100"
        assert rows[1][2] == "1000" and rows[1][12] == "80"
        # Second position, highest score first
        assert rows[2][2] == "2000" and rows[2][12] == "90"
        assert rows[3][2] == "2000" and rows[3][12] == "50"

    def test_empty_hits(self, tmp_path):
        output = str(tmp_path / "test.tsv")
        write_tsv([], output)

        with open(output) as fh:
            lines = fh.readlines()
        # Only header
        assert len(lines) == 1
