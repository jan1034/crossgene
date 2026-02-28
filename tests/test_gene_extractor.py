"""Tests for gene_extractor module."""

import logging
import textwrap
from pathlib import Path

import pysam
import pytest

from compare_genes.gene_extractor import extract_sequence, lookup_gene, _reverse_complement

TEST_DATA = Path(__file__).parent / "test_data"


@pytest.fixture
def tiny_gtf(tmp_path):
    """Create a minimal GTF with two genes (one plus, one minus strand)."""
    gtf = tmp_path / "test.gtf"
    # GTF coordinates are 1-based inclusive
    gtf.write_text(textwrap.dedent("""\
        chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G001"; gene_name "GENE_A"; gene_biotype "protein_coding";
        chr1\ttest\tgene\t301\t400\t.\t-\t.\tgene_id "G002"; gene_name "GENE_B"; gene_biotype "protein_coding";
        chr2\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G003"; gene_name "DUP_GENE"; gene_biotype "protein_coding";
        chr3\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G004"; gene_name "DUP_GENE"; gene_biotype "protein_coding";
    """))
    return str(gtf)


@pytest.fixture
def tiny_fasta(tmp_path):
    """Create a minimal FASTA with known sequences."""
    fasta_path = tmp_path / "test.fa"
    # chr1: 500 bases, chr2: 300 bases, chr3: 300 bases
    # GENE_A is at chr1:101-200 (1-based) → 0-based [100, 200) = 100 bases
    # GENE_B is at chr1:301-400 (1-based) → 0-based [300, 400) = 100 bases, minus strand
    chr1_seq = "N" * 100 + "ACGT" * 25 + "N" * 100 + "AAAA" * 25 + "N" * 100
    chr2_seq = "N" * 100 + "CCCC" * 25 + "N" * 100
    chr3_seq = "N" * 100 + "GGGG" * 25 + "N" * 100
    fasta_path.write_text(f">chr1\n{chr1_seq}\n>chr2\n{chr2_seq}\n>chr3\n{chr3_seq}\n")
    # Create index
    pysam.faidx(str(fasta_path))
    return str(fasta_path)


class TestLookupGene:
    def test_finds_gene(self, tiny_gtf):
        gene = lookup_gene("GENE_A", tiny_gtf)
        assert gene.name == "GENE_A"
        assert gene.gene_id == "G001"
        assert gene.chrom == "chr1"
        assert gene.start == 100  # 0-based
        assert gene.end == 200
        assert gene.strand == "+"
        assert gene.sequence == ""

    def test_minus_strand_gene(self, tiny_gtf):
        gene = lookup_gene("GENE_B", tiny_gtf)
        assert gene.strand == "-"
        assert gene.start == 300  # 0-based
        assert gene.end == 400

    def test_gene_not_found(self, tiny_gtf):
        with pytest.raises(ValueError, match="not found"):
            lookup_gene("NONEXISTENT", tiny_gtf)

    def test_multiple_matches_warns(self, tiny_gtf, caplog):
        with caplog.at_level(logging.WARNING):
            gene = lookup_gene("DUP_GENE", tiny_gtf)
        assert "Multiple matches" in caplog.text
        # Should return the first match
        assert gene.gene_id == "G003"
        assert gene.chrom == "chr2"


class TestExtractSequence:
    def test_plus_strand(self, tiny_gtf, tiny_fasta):
        gene = lookup_gene("GENE_A", tiny_gtf)
        gene = extract_sequence(gene, tiny_fasta)
        assert len(gene.sequence) == 100
        assert gene.sequence == "ACGT" * 25

    def test_minus_strand_revcomp(self, tiny_gtf, tiny_fasta):
        gene = lookup_gene("GENE_B", tiny_gtf)
        gene = extract_sequence(gene, tiny_fasta)
        assert len(gene.sequence) == 100
        # Original: AAAA * 25 = 100 A's
        # Reverse complement of all A's = all T's
        assert gene.sequence == "TTTT" * 25


class TestReverseComplement:
    def test_basic(self):
        assert _reverse_complement("ACGT") == "ACGT"  # palindromic
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("ATCG") == "CGAT"

    def test_lowercase(self):
        assert _reverse_complement("acgt") == "acgt"
