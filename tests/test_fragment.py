"""Tests for fragment module."""

from pathlib import Path

import pytest

from crossgene.fragment import fragment_index_to_genomic, generate_fragments
from crossgene.models import GeneRecord


def _make_gene(seq: str, strand: str = "+", start: int = 1000, name: str = "TEST") -> GeneRecord:
    """Helper to create a GeneRecord with a given sequence."""
    return GeneRecord(
        name=name,
        gene_id="G001",
        chrom="chr1",
        start=start,
        end=start + len(seq),
        strand=strand,
        sequence=seq,
    )


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA file and return list of (header, sequence) tuples."""
    entries = []
    header = None
    seq_lines = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_lines)))
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


class TestGenerateFragments:
    def test_basic_fragmentation(self):
        """500bp sequence, fragment_size=50, step=10 → known fragment count."""
        seq = "A" * 500
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=50, step_size=10)
        entries = _read_fasta(path)

        # Full fragments: positions 0,10,20,...,450 → (500-50)//10 + 1 = 46
        # Position 460: remainder = 40bp < 50/2=25? No, 40 >= 25 → included
        assert len(entries) == 46 + 1  # 46 full + 1 partial
        path.unlink()

    def test_fragment_naming(self):
        seq = "ACGT" * 25  # 100bp
        gene = _make_gene(seq, name="BRCA1")
        path = generate_fragments(gene, fragment_size=50, step_size=50)
        entries = _read_fasta(path)

        assert entries[0][0] == "BRCA1:0"
        assert entries[1][0] == "BRCA1:1"
        assert len(entries) == 2
        path.unlink()

    def test_fragment_content(self):
        seq = "AAAA" + "CCCC" + "GGGG"  # 12bp
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=4, step_size=4)
        entries = _read_fasta(path)

        assert entries[0][1] == "AAAA"
        assert entries[1][1] == "CCCC"
        assert entries[2][1] == "GGGG"
        assert len(entries) == 3
        path.unlink()

    def test_short_sequence(self):
        """Sequence shorter than fragment_size → single fragment."""
        seq = "ACGT"
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=50, step_size=1)
        entries = _read_fasta(path)

        assert len(entries) == 1
        assert entries[0][1] == "ACGT"
        path.unlink()

    def test_trailing_fragment_included(self):
        """Trailing fragment >= fragment_size/2 should be included."""
        seq = "A" * 75  # frag=50, step=50 → full at 0, remainder at 50 = 25bp = 50/2 → include
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=50, step_size=50)
        entries = _read_fasta(path)

        assert len(entries) == 2
        assert len(entries[1][1]) == 25
        path.unlink()

    def test_trailing_fragment_skipped(self):
        """Trailing fragment < fragment_size/2 should be skipped."""
        seq = "A" * 74  # frag=50, step=50 → full at 0, remainder at 50 = 24bp < 25 → skip
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=50, step_size=50)
        entries = _read_fasta(path)

        assert len(entries) == 1
        path.unlink()

    def test_step_greater_than_one(self):
        seq = "A" * 100
        gene = _make_gene(seq)
        path = generate_fragments(gene, fragment_size=50, step_size=25)
        entries = _read_fasta(path)

        # Full: pos 0, 25, 50 → 3 fragments
        # Partial: pos 75 → 25bp = fragment_size/2 → included
        assert len(entries) == 4
        path.unlink()

    def test_empty_sequence_raises(self):
        gene = _make_gene("")
        with pytest.raises(ValueError, match="no sequence"):
            generate_fragments(gene, fragment_size=50, step_size=1)


class TestFragmentIndexToGenomic:
    def test_plus_strand(self):
        gene = _make_gene("A" * 100, strand="+", start=1000)
        # Fragment 0: sense [0,50) → genomic [1000, 1050)
        assert fragment_index_to_genomic(gene, 0, 50, 10) == (1000, 1050)
        # Fragment 1: sense [10,60) → genomic [1010, 1060)
        assert fragment_index_to_genomic(gene, 1, 50, 10) == (1010, 1060)

    def test_minus_strand(self):
        gene = _make_gene("A" * 100, strand="-", start=1000)
        # gene.end = 1100
        # Fragment 0: sense [0,50) → genomic [1100-0-50, 1100-0) = [1050, 1100)
        assert fragment_index_to_genomic(gene, 0, 50, 10) == (1050, 1100)
        # Fragment 1: sense [10,60) → genomic [1100-10-50, 1100-10) = [1040, 1090)
        assert fragment_index_to_genomic(gene, 1, 50, 10) == (1040, 1090)

    def test_trailing_partial_fragment(self):
        gene = _make_gene("A" * 75, strand="+", start=1000)
        # Fragment at sense pos 50: only 25bp remain
        # idx=1 with step=50: sense [50,75) → genomic [1050, 1075)
        assert fragment_index_to_genomic(gene, 1, 50, 50) == (1050, 1075)

    def test_trailing_partial_minus_strand(self):
        gene = _make_gene("A" * 75, strand="-", start=1000)
        # gene.end = 1075
        # idx=1, step=50: sense [50,75) → 25bp
        # genomic: [1075-50-25, 1075-50) = [1000, 1025)
        assert fragment_index_to_genomic(gene, 1, 50, 50) == (1000, 1025)
