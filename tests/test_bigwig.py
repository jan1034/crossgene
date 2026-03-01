"""Tests for bigwig module."""

import numpy as np
import pyBigWig
import pytest

from crossgene.bigwig import read_chrom_sizes, write_bigwig


@pytest.fixture
def chrom_sizes_file(tmp_path):
    """Create a minimal chrom.sizes file."""
    path = tmp_path / "chrom.sizes"
    path.write_text("chr1\t1000000\nchr2\t800000\n")
    return str(path)


@pytest.fixture
def chrom_sizes_with_header(tmp_path):
    """Create a chrom.sizes file with a header line (like the reference file)."""
    path = tmp_path / "chrom.sizes"
    path.write_text("id\tlength\tlength\nchr1\t1000000\t1000000\nchr2\t800000\t800000\n")
    return str(path)


class TestReadChromSizes:
    def test_basic(self, chrom_sizes_file):
        sizes = read_chrom_sizes(chrom_sizes_file)
        assert sizes == {"chr1": 1000000, "chr2": 800000}

    def test_with_header(self, chrom_sizes_with_header):
        sizes = read_chrom_sizes(chrom_sizes_with_header)
        assert sizes == {"chr1": 1000000, "chr2": 800000}


class TestWriteBigwig:
    def test_write_and_read_back(self, tmp_path, chrom_sizes_file):
        """Write scores, read back with pyBigWig, verify values."""
        scores = np.array([100.0, 95.0, 80.0, 0.0, 50.0], dtype=np.float64)
        output = str(tmp_path / "test.bw")

        write_bigwig(scores, "chr1", 5000, chrom_sizes_file, output)

        bw = pyBigWig.open(output)
        # Read values back at each position
        for i, expected in enumerate(scores):
            vals = bw.values("chr1", 5000 + i, 5000 + i + 1)
            assert len(vals) == 1
            assert abs(vals[0] - expected) < 0.01, f"Position {5000+i}: {vals[0]} != {expected}"
        bw.close()

    def test_larger_array(self, tmp_path, chrom_sizes_file):
        """Write a 1000-element score array."""
        scores = np.random.default_rng(42).uniform(0, 100, size=1000)
        output = str(tmp_path / "test.bw")

        write_bigwig(scores, "chr1", 10000, chrom_sizes_file, output)

        bw = pyBigWig.open(output)
        vals = bw.values("chr1", 10000, 11000)
        np.testing.assert_allclose(vals, scores, atol=0.01)
        bw.close()

    def test_unknown_chrom_raises(self, tmp_path, chrom_sizes_file):
        scores = np.array([50.0])
        output = str(tmp_path / "test.bw")
        with pytest.raises(ValueError, match="not found"):
            write_bigwig(scores, "chrZ", 0, chrom_sizes_file, output)

    def test_empty_scores(self, tmp_path, chrom_sizes_file):
        """Empty score array produces a valid (empty) BigWig."""
        scores = np.array([], dtype=np.float64)
        output = str(tmp_path / "test.bw")
        write_bigwig(scores, "chr1", 0, chrom_sizes_file, output)

        bw = pyBigWig.open(output)
        # No entries, but file is valid
        assert bw.chroms() is not None
        bw.close()
