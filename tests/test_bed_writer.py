"""Tests for the BED9 hit export writer."""

from __future__ import annotations

import os

import pytest

from crossgene.bed_writer import COLOR_ANTISENSE, COLOR_SENSE, write_hit_beds
from crossgene.models import AlignmentHit, GeneRecord


def _make_gene(name: str = "GENEA", chrom: str = "chr1", start: int = 1000, end: int = 2000) -> GeneRecord:
    return GeneRecord(
        name=name, gene_id=f"ENSG_{name}", chrom=chrom,
        start=start, end=end, strand="+", sequence="A" * (end - start),
    )


def _make_hit(
    query_start: int = 100, query_end: int = 150,
    target_start: int = 500, target_end: int = 550,
    strand: str = "+", identity: float = 0.90, is_primary: bool = True,
) -> AlignmentHit:
    return AlignmentHit(
        query_chrom="chr1", query_start=query_start, query_end=query_end,
        target_chrom="chr2", target_start=target_start, target_end=target_end,
        strand=strand, identity=identity, query_coverage=1.0, mapq=60,
        cigar="50M", alignment_score=100, is_primary=is_primary,
        query_gene="GENEA", target_gene="GENEB", direction="A→B",
    )


def _parse_bed(path: str) -> tuple[str, list[list[str]]]:
    with open(path) as f:
        lines = f.read().strip().split("\n")
    header = lines[0] if lines else ""
    rows = [line.split("\t") for line in lines[1:]] if len(lines) > 1 else []
    return header, rows


class TestNumberingAndCrossrefs:
    def test_numbering_and_crossrefs(self, tmp_path):
        hits = [
            _make_hit(query_start=300, query_end=350, target_start=800, target_end=850, identity=0.9),
            _make_hit(query_start=100, query_end=150, target_start=600, target_end=650, identity=0.95),
            _make_hit(query_start=200, query_end=250, target_start=500, target_end=550, identity=0.85),
        ]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB", chrom="chr2"), "AtoB", str(tmp_path))
        _, q_rows = _parse_bed(q_bed)
        _, t_rows = _parse_bed(t_bed)

        assert [r[1] for r in q_rows] == ["100", "200", "300"]
        q_names = [r[3] for r in q_rows]
        assert q_names[0].startswith("A1-") and q_names[2].startswith("A3-")

        q_crossrefs = {n.split("-")[0]: n.split("-")[1] for n in q_names}
        t_crossrefs = {n.split("-")[0]: n.split("-")[1] for n in [r[3] for r in t_rows]}
        for qk, tv in q_crossrefs.items():
            assert t_crossrefs[tv] == qk


class TestBed9Format:
    def test_format_columns_score_thick_rgb(self, tmp_path):
        """Verify 9 columns, score 0-1000, thickStart/End == start/end, and strand coloring."""
        hits = [_make_hit(strand="+", identity=0.85), _make_hit(query_start=200, strand="-", identity=1.0)]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))

        for bed_path in (q_bed, t_bed):
            _, rows = _parse_bed(bed_path)
            for row in rows:
                assert len(row) == 9
                assert 0 <= int(row[4]) <= 1000
                assert row[6] == row[1]  # thickStart == start
                assert row[7] == row[2]  # thickEnd == end

        _, q_rows = _parse_bed(q_bed)
        rgbs = {row[8] for row in q_rows}
        assert COLOR_SENSE in rgbs and COLOR_ANTISENSE in rgbs


class TestTrackHeaderAndFilters:
    def test_track_header(self, tmp_path):
        hits = [_make_hit()]
        q_bed, _ = write_hit_beds(hits, _make_gene("BRCA1"), _make_gene("BRCA2"), "AtoB", str(tmp_path))
        header, _ = _parse_bed(q_bed)
        assert header.startswith("track name=")
        assert "BRCA1" in header

    def test_primary_only_vs_all_hits(self, tmp_path):
        hits = [_make_hit(is_primary=True), _make_hit(query_start=200, is_primary=False)]

        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        assert len(rows) == 1  # default: primary only

        q_bed2, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path), primary_only=False)
        _, rows2 = _parse_bed(q_bed2)
        assert len(rows2) == 2


class TestEdgeCases:
    def test_empty_hits_no_crash(self, tmp_path):
        q_bed, t_bed = write_hit_beds([], _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        assert os.path.exists(q_bed) and os.path.exists(t_bed)
        _, rows = _parse_bed(q_bed)
        assert len(rows) == 0

    @pytest.mark.parametrize("dir_tag,q_prefix,t_prefix", [
        ("AtoB", "A", "B"),
        ("BtoA", "B", "A"),
    ])
    def test_direction_prefixes(self, tmp_path, dir_tag, q_prefix, t_prefix):
        hits = [_make_hit()]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), dir_tag, str(tmp_path))
        _, q_rows = _parse_bed(q_bed)
        _, t_rows = _parse_bed(t_bed)
        assert q_rows[0][3].startswith(q_prefix)
        assert t_rows[0][3].startswith(t_prefix)
