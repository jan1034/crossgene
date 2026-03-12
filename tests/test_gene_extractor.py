"""Tests for gene_extractor module."""

import logging
import textwrap
from pathlib import Path

import pysam
import pytest

from crossgene.gene_extractor import extract_sequence, load_features, lookup_gene
from crossgene.models import GeneRecord

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

    def test_minus_strand_no_revcomp(self, tiny_gtf, tiny_fasta):
        gene = lookup_gene("GENE_B", tiny_gtf)
        gene = extract_sequence(gene, tiny_fasta)
        assert len(gene.sequence) == 100
        # Genomic orientation: raw sequence is all A's
        assert gene.sequence == "AAAA" * 25


@pytest.fixture
def annotation_gtf(tmp_path):
    """Create an annotation GTF with transcripts, exons, and CDS for GENE_A (G001).

    Two transcripts:
      - T001: tagged Ensembl_canonical, has 2 exons + 1 CDS
      - T002: non-canonical, has 1 exon + 1 CDS
    """
    gtf = tmp_path / "annotation.gtf"
    # GTF coordinates are 1-based inclusive
    gtf.write_text(textwrap.dedent("""\
        chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G001"; gene_name "GENE_A";
        chr1\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "G001"; transcript_id "T001"; tag "basic,Ensembl_canonical";
        chr1\ttest\texon\t101\t130\t.\t+\t.\tgene_id "G001"; transcript_id "T001";
        chr1\ttest\texon\t161\t200\t.\t+\t.\tgene_id "G001"; transcript_id "T001";
        chr1\ttest\tCDS\t105\t125\t.\t+\t.\tgene_id "G001"; transcript_id "T001";
        chr1\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "G001"; transcript_id "T002"; tag "basic";
        chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "G001"; transcript_id "T002";
        chr1\ttest\tCDS\t110\t145\t.\t+\t.\tgene_id "G001"; transcript_id "T002";
    """))
    return str(gtf)


def _make_gene_record(name="GENE_A", gene_id="G001", chrom="chr1",
                       start=100, end=200, strand="+"):
    """Helper to create a bare GeneRecord for load_features tests."""
    return GeneRecord(
        name=name, gene_id=gene_id, chrom=chrom,
        start=start, end=end, strand=strand, sequence="",
    )


class TestLoadFeatures:
    def test_canonical_mode(self, annotation_gtf):
        """Canonical mode should return only T001 features (2 exons + 1 CDS)."""
        gene = _make_gene_record()
        result = load_features(gene, annotation_gtf, transcript_mode="canonical")
        assert len(result.features) == 3
        types = [f.feature_type for f in result.features]
        assert types.count("exon") == 2
        assert types.count("CDS") == 1
        # First exon: GTF 101-130 → 0-based [100, 130)
        assert result.features[0].start == 100
        assert result.features[0].end == 130

    def test_all_mode(self, annotation_gtf):
        """All mode should return features from both transcripts, deduplicated."""
        gene = _make_gene_record()
        result = load_features(gene, annotation_gtf, transcript_mode="all")
        # T001: exon(101,130), exon(161,200), CDS(105,125)
        # T002: exon(101,150), CDS(110,145)
        # All have unique (type, start, end) so 5 total
        assert len(result.features) == 5

    def test_custom_feature_types(self, annotation_gtf):
        """Only requested feature types should be returned."""
        gene = _make_gene_record()
        result = load_features(
            gene, annotation_gtf, transcript_mode="canonical",
            feature_types={"exon"},
        )
        assert all(f.feature_type == "exon" for f in result.features)
        assert len(result.features) == 2

    def test_gene_not_in_annotation(self, annotation_gtf):
        """Missing gene_id should return original record unchanged."""
        gene = _make_gene_record(gene_id="MISSING")
        result = load_features(gene, annotation_gtf)
        assert result is gene
        assert result.features == []

    def test_preserves_existing_fields(self, annotation_gtf):
        """load_features should preserve name, chrom, strand, sequence, etc."""
        gene = _make_gene_record()
        gene = GeneRecord(
            name=gene.name, gene_id=gene.gene_id, chrom=gene.chrom,
            start=gene.start, end=gene.end, strand=gene.strand,
            sequence="ACGT" * 25,
        )
        result = load_features(gene, annotation_gtf)
        assert result.sequence == "ACGT" * 25
        assert result.name == "GENE_A"
        assert result.chrom == "chr1"
        assert len(result.features) > 0

    def test_features_sorted_by_start(self, annotation_gtf):
        """Features should be sorted by start coordinate."""
        gene = _make_gene_record()
        result = load_features(gene, annotation_gtf, transcript_mode="canonical")
        starts = [f.start for f in result.features]
        assert starts == sorted(starts)


