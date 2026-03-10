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
    query_start: int = 100,
    query_end: int = 150,
    target_start: int = 500,
    target_end: int = 550,
    strand: str = "+",
    identity: float = 0.90,
    is_primary: bool = True,
    query_chrom: str = "chr1",
    target_chrom: str = "chr2",
) -> AlignmentHit:
    return AlignmentHit(
        query_chrom=query_chrom, query_start=query_start, query_end=query_end,
        target_chrom=target_chrom, target_start=target_start, target_end=target_end,
        strand=strand, identity=identity, query_coverage=1.0, mapq=60,
        cigar="50M", alignment_score=100, is_primary=is_primary,
        query_gene="GENEA", target_gene="GENEB", direction="A→B",
    )


def _parse_bed(path: str) -> tuple[str, list[list[str]]]:
    """Return (header_line, list_of_row_columns) from a BED file."""
    with open(path) as f:
        lines = f.read().strip().split("\n")
    header = lines[0] if lines else ""
    rows = [line.split("\t") for line in lines[1:]] if len(lines) > 1 else []
    return header, rows


class TestNumberingAndCrossrefs:
    def test_numbering_and_crossrefs(self, tmp_path):
        """Query-side and target-side numbers are assigned by position, cross-refs match."""
        hits = [
            _make_hit(query_start=300, query_end=350, target_start=800, target_end=850, identity=0.9),
            _make_hit(query_start=100, query_end=150, target_start=600, target_end=650, identity=0.95),
            _make_hit(query_start=200, query_end=250, target_start=500, target_end=550, identity=0.85),
        ]
        q_gene = _make_gene("GENEA")
        t_gene = _make_gene("GENEB", chrom="chr2")

        q_bed, t_bed = write_hit_beds(hits, q_gene, t_gene, "AtoB", str(tmp_path))

        _, q_rows = _parse_bed(q_bed)
        _, t_rows = _parse_bed(t_bed)

        # Query BED sorted by query_start: 100, 200, 300 → A1, A2, A3
        assert [r[1] for r in q_rows] == ["100", "200", "300"]
        q_names = [r[3] for r in q_rows]
        assert q_names[0].startswith("A1-")
        assert q_names[1].startswith("A2-")
        assert q_names[2].startswith("A3-")

        # Target BED sorted by target_start: 500, 600, 800 → B1, B2, B3
        assert [r[1] for r in t_rows] == ["500", "600", "800"]
        t_names = [r[3] for r in t_rows]
        assert t_names[0].startswith("B1-")
        assert t_names[1].startswith("B2-")
        assert t_names[2].startswith("B3-")

        # Cross-references match: if query has A1-Bx, target must have Bx-A1
        q_crossrefs = {n.split("-")[0]: n.split("-")[1] for n in q_names}
        t_crossrefs = {n.split("-")[0]: n.split("-")[1] for n in t_names}
        for qk, tv in q_crossrefs.items():
            assert t_crossrefs[tv] == qk


class TestBed9Format:
    def test_bed9_columns(self, tmp_path):
        """Each data row has exactly 9 tab-separated columns."""
        hits = [_make_hit(identity=0.85), _make_hit(query_start=200, identity=1.0)]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))

        for bed_path in (q_bed, t_bed):
            _, rows = _parse_bed(bed_path)
            for row in rows:
                assert len(row) == 9

    def test_score_range(self, tmp_path):
        hits = [_make_hit(identity=0.0), _make_hit(query_start=200, identity=1.0)]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        for row in rows:
            score = int(row[4])
            assert 0 <= score <= 1000

    def test_thick_equals_coords(self, tmp_path):
        hits = [_make_hit(query_start=100, query_end=150)]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, q_rows = _parse_bed(q_bed)
        for row in q_rows:
            assert row[6] == row[1]  # thickStart == start
            assert row[7] == row[2]  # thickEnd == end

    def test_itemrgb_values(self, tmp_path):
        hits = [_make_hit(strand="+"), _make_hit(query_start=200, strand="-")]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        rgbs = {row[8] for row in rows}
        assert rgbs == {COLOR_SENSE, COLOR_ANTISENSE}


class TestTrackHeader:
    def test_track_header(self, tmp_path):
        hits = [_make_hit()]
        q_bed, t_bed = write_hit_beds(hits, _make_gene("BRCA1"), _make_gene("BRCA2"), "AtoB", str(tmp_path))
        header, _ = _parse_bed(q_bed)
        assert header.startswith("track name=")
        assert "BRCA1" in header
        assert "AtoB" in header


class TestPrimaryOnlyFilter:
    def test_default_excludes_secondary(self, tmp_path):
        hits = [
            _make_hit(is_primary=True),
            _make_hit(query_start=200, is_primary=False),
        ]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        assert len(rows) == 1

    def test_all_hits_includes_secondary(self, tmp_path):
        hits = [
            _make_hit(is_primary=True),
            _make_hit(query_start=200, is_primary=False),
        ]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path), primary_only=False)
        _, rows = _parse_bed(q_bed)
        assert len(rows) == 2


class TestStrandColoring:
    def test_sense_steelblue(self, tmp_path):
        hits = [_make_hit(strand="+")]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        assert rows[0][8] == COLOR_SENSE

    def test_antisense_firebrick(self, tmp_path):
        hits = [_make_hit(strand="-")]
        q_bed, _ = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, rows = _parse_bed(q_bed)
        assert rows[0][8] == COLOR_ANTISENSE


class TestEmptyHits:
    def test_empty_hits_no_crash(self, tmp_path):
        q_bed, t_bed = write_hit_beds([], _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        assert os.path.exists(q_bed)
        assert os.path.exists(t_bed)
        header, rows = _parse_bed(q_bed)
        assert header.startswith("track name=")
        assert len(rows) == 0


class TestDirectionPrefixes:
    def test_atob_prefixes(self, tmp_path):
        hits = [_make_hit()]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "AtoB", str(tmp_path))
        _, q_rows = _parse_bed(q_bed)
        _, t_rows = _parse_bed(t_bed)
        assert q_rows[0][3].startswith("A")
        assert t_rows[0][3].startswith("B")

    def test_btoa_prefixes(self, tmp_path):
        hits = [_make_hit()]
        q_bed, t_bed = write_hit_beds(hits, _make_gene(), _make_gene("GENEB"), "BtoA", str(tmp_path))
        _, q_rows = _parse_bed(q_bed)
        _, t_rows = _parse_bed(t_bed)
        assert q_rows[0][3].startswith("B")
        assert t_rows[0][3].startswith("A")
