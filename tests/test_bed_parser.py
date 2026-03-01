"""Tests for bed_parser module."""

import logging

import pytest

from crossgene.bed_parser import filter_and_clip, parse_bed
from crossgene.models import BedRegion


class TestParseBed:
    def test_bed6(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\tAluSc\t2224\t-\nchr2\t300\t400\tL1\t500\t+\n")
        regions = parse_bed(str(bed))
        assert len(regions) == 2
        assert regions[0] == BedRegion("chr1", 100, 200, "AluSc", 2224, "-")
        assert regions[1] == BedRegion("chr2", 300, 400, "L1", 500, "+")

    def test_bed3_defaults(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\n")
        regions = parse_bed(str(bed))
        assert len(regions) == 1
        assert regions[0].name == "."
        assert regions[0].score == 0
        assert regions[0].strand == "."

    def test_bed4(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\tMyRegion\n")
        regions = parse_bed(str(bed))
        assert regions[0].name == "MyRegion"
        assert regions[0].score == 0

    def test_skip_comments_track_browser_empty(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("# comment\ntrack name=test\nbrowser position\n\nchr1\t100\t200\n")
        regions = parse_bed(str(bed))
        assert len(regions) == 1

    def test_malformed_too_few_columns(self, tmp_path, caplog):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\n")
        with caplog.at_level(logging.WARNING):
            regions = parse_bed(str(bed))
        assert len(regions) == 0
        assert "at least 3 columns" in caplog.text

    def test_malformed_non_numeric(self, tmp_path, caplog):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\tabc\t200\n")
        with caplog.at_level(logging.WARNING):
            regions = parse_bed(str(bed))
        assert len(regions) == 0
        assert "non-numeric" in caplog.text

    def test_malformed_start_ge_end(self, tmp_path, caplog):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t200\t200\tEqual\n")
        with caplog.at_level(logging.WARNING):
            regions = parse_bed(str(bed))
        assert len(regions) == 0
        assert "start >= end" in caplog.text

    def test_non_numeric_score_defaults(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\tName\tnotanumber\t+\n")
        regions = parse_bed(str(bed))
        assert regions[0].score == 0


class TestFilterAndClip:
    def _regions(self):
        return [
            BedRegion("chr1", 100, 200, "A"),
            BedRegion("chr1", 180, 300, "B"),
            BedRegion("chr1", 400, 500, "C"),
            BedRegion("chr2", 100, 200, "D"),
        ]

    def test_fully_inside(self):
        regions = [BedRegion("chr1", 120, 180, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 1
        assert result[0].start == 120
        assert result[0].end == 180

    def test_partially_overlapping_left(self):
        regions = [BedRegion("chr1", 50, 150, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert result[0].start == 100
        assert result[0].end == 150

    def test_partially_overlapping_right(self):
        regions = [BedRegion("chr1", 150, 250, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert result[0].start == 150
        assert result[0].end == 200

    def test_fully_outside(self):
        regions = [BedRegion("chr1", 300, 400, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 0

    def test_wrong_chromosome(self):
        regions = [BedRegion("chr2", 100, 200, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 0

    def test_edge_abutting_no_overlap(self):
        """Region [200,300) does not overlap [100,200)."""
        regions = [BedRegion("chr1", 200, 300, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 0

    def test_does_not_mutate_input(self):
        original = BedRegion("chr1", 50, 250, "X")
        regions = [original]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert original.start == 50  # unchanged
        assert original.end == 250
        assert result[0].start == 100
        assert result[0].end == 200

    def test_multiple_regions_mixed(self):
        result = filter_and_clip(self._regions(), "chr1", 150, 350)
        names = [r.name for r in result]
        assert "A" in names  # clipped to 150-200
        assert "B" in names  # clipped to 180-300
        assert "C" not in names  # 400-500 outside
        assert "D" not in names  # wrong chrom
