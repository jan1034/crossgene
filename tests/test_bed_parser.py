"""Tests for bed_parser module."""

import logging

import pytest

from crossgene.bed_parser import filter_and_clip, parse_bed
from crossgene.models import BedRegion


class TestParseBed:
    @pytest.mark.parametrize("content,expected_name,expected_score,expected_strand", [
        ("chr1\t100\t200\n", ".", 0, "."),                              # BED3 defaults
        ("chr1\t100\t200\tMyRegion\n", "MyRegion", 0, "."),             # BED4
        ("chr1\t100\t200\tAluSc\t2224\t-\n", "AluSc", 2224, "-"),      # BED6
    ])
    def test_bed_formats(self, tmp_path, content, expected_name, expected_score, expected_strand):
        bed = tmp_path / "test.bed"
        bed.write_text(content)
        regions = parse_bed(str(bed))
        assert len(regions) == 1
        assert regions[0].chrom == "chr1"
        assert regions[0].start == 100
        assert regions[0].end == 200
        assert regions[0].name == expected_name
        assert regions[0].score == expected_score
        assert regions[0].strand == expected_strand

    def test_skip_comments_track_browser_empty(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("# comment\ntrack name=test\nbrowser position\n\nchr1\t100\t200\n")
        regions = parse_bed(str(bed))
        assert len(regions) == 1

    @pytest.mark.parametrize("content,expected_warning", [
        ("chr1\t100\n", "at least 3 columns"),
        ("chr1\tabc\t200\n", "non-numeric"),
        ("chr1\t200\t200\tEqual\n", "start >= end"),
    ])
    def test_malformed_lines(self, tmp_path, caplog, content, expected_warning):
        bed = tmp_path / "test.bed"
        bed.write_text(content)
        with caplog.at_level(logging.WARNING):
            regions = parse_bed(str(bed))
        assert len(regions) == 0
        assert expected_warning in caplog.text

    def test_non_numeric_score_defaults(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\tName\tnotanumber\t+\n")
        regions = parse_bed(str(bed))
        assert regions[0].score == 0


class TestFilterAndClip:
    def test_clipping(self):
        """Tests fully inside, left overlap, and right overlap clipping."""
        regions = [
            BedRegion("chr1", 120, 180, "inside"),      # fully inside
            BedRegion("chr1", 50, 150, "left"),          # left overlap
            BedRegion("chr1", 150, 250, "right"),        # right overlap
        ]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 3
        assert (result[0].start, result[0].end) == (120, 180)  # unchanged
        assert (result[1].start, result[1].end) == (100, 150)  # clipped left
        assert (result[2].start, result[2].end) == (150, 200)  # clipped right

    @pytest.mark.parametrize("chrom,start,end", [
        ("chr1", 300, 400),   # fully outside
        ("chr2", 100, 200),   # wrong chromosome
        ("chr1", 200, 300),   # abutting, no overlap
    ])
    def test_no_overlap(self, chrom, start, end):
        regions = [BedRegion(chrom, start, end, "X")]
        result = filter_and_clip(regions, "chr1", 100, 200)
        assert len(result) == 0

    def test_does_not_mutate_input(self):
        original = BedRegion("chr1", 50, 250, "X")
        result = filter_and_clip([original], "chr1", 100, 200)
        assert original.start == 50  # unchanged
        assert result[0].start == 100

    def test_multiple_regions_mixed(self):
        regions = [
            BedRegion("chr1", 100, 200, "A"),
            BedRegion("chr1", 180, 300, "B"),
            BedRegion("chr1", 400, 500, "C"),
            BedRegion("chr2", 100, 200, "D"),
        ]
        result = filter_and_clip(regions, "chr1", 150, 350)
        names = [r.name for r in result]
        assert "A" in names and "B" in names
        assert "C" not in names and "D" not in names
