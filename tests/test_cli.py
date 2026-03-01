"""Tests for CLI entry point."""

import textwrap

import pysam
import pytest
from click.testing import CliRunner

from crossgene.cli import main, _parse_formats


class TestParseFormats:
    def test_all_formats(self):
        assert _parse_formats("bigwig,tsv,plot") == {"bigwig", "tsv", "plot"}

    def test_single_format(self):
        assert _parse_formats("tsv") == {"tsv"}

    def test_invalid_format(self):
        with pytest.raises(Exception):
            _parse_formats("bigwig,invalid")


class TestCLI:
    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "Compare two gene sequences" in result.output

    def test_missing_required(self):
        runner = CliRunner()
        result = runner.invoke(main, [])
        assert result.exit_code != 0

    def test_gene_not_found(self, tmp_path):
        """Gene not found gives a clear error, not a traceback."""
        gtf = tmp_path / "test.gtf"
        gtf.write_text(textwrap.dedent("""\
            chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G001"; gene_name "REAL_GENE"; gene_biotype "protein_coding";
        """))
        genome = tmp_path / "test.fa"
        genome.write_text(">chr1\n" + "A" * 300 + "\n")
        pysam.faidx(str(genome))
        chrom_sizes = tmp_path / "chrom.sizes"
        chrom_sizes.write_text("chr1\t300\n")

        runner = CliRunner()
        result = runner.invoke(main, [
            "--gene-a", "NONEXISTENT",
            "--gene-b", "REAL_GENE",
            "--gtf", str(gtf),
            "--genome", str(genome),
            "--chrom-sizes", str(chrom_sizes),
            "--outdir", str(tmp_path / "out"),
        ])
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_integration_synthetic(self, tmp_path):
        """End-to-end with tiny synthetic genes."""
        # Create synthetic GTF with two genes
        shared = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 52bp
        seq_a = "N" * 50 + shared + "N" * 48  # 150bp total
        seq_b = "T" * 30 + shared + "T" * 68  # 150bp total

        gtf = tmp_path / "test.gtf"
        gtf.write_text(textwrap.dedent("""\
            chr1\ttest\tgene\t1\t150\t.\t+\t.\tgene_id "G001"; gene_name "SYNTH_A"; gene_biotype "protein_coding";
            chr2\ttest\tgene\t1\t150\t.\t+\t.\tgene_id "G002"; gene_name "SYNTH_B"; gene_biotype "protein_coding";
        """))

        genome = tmp_path / "test.fa"
        genome.write_text(f">chr1\n{seq_a}\n>chr2\n{seq_b}\n")
        pysam.faidx(str(genome))

        chrom_sizes = tmp_path / "chrom.sizes"
        chrom_sizes.write_text("chr1\t150\nchr2\t150\n")

        outdir = tmp_path / "output"
        runner = CliRunner()
        result = runner.invoke(main, [
            "--gene-a", "SYNTH_A",
            "--gene-b", "SYNTH_B",
            "--fragment-size", "50",
            "--step-size", "25",
            "--min-quality", "30",
            "--gtf", str(gtf),
            "--genome", str(genome),
            "--chrom-sizes", str(chrom_sizes),
            "--outdir", str(outdir),
            "--output-formats", "bigwig,tsv,plot",
        ])

        assert result.exit_code == 0, f"CLI failed:\n{result.output}"

        # Check all 5 output files exist
        assert (outdir / "SYNTH_A_vs_SYNTH_B.mappability.bw").exists()
        assert (outdir / "SYNTH_B_vs_SYNTH_A.mappability.bw").exists()
        assert (outdir / "SYNTH_A_vs_SYNTH_B.hits.tsv").exists()
        assert (outdir / "SYNTH_B_vs_SYNTH_A.hits.tsv").exists()
        assert (outdir / "SYNTH_A_vs_SYNTH_B.circlize.pdf").exists()

        # TSV should have a header + at least some hits
        tsv_lines = (outdir / "SYNTH_A_vs_SYNTH_B.hits.tsv").read_text().strip().splitlines()
        assert len(tsv_lines) >= 2  # header + at least 1 hit
