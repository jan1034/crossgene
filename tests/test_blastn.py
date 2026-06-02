"""Tests for blastn module."""

from unittest.mock import MagicMock, patch

import pytest

from crossgene.blastn import (
    BlastParams,
    _build_blastn_command,
    _make_blast_db,
    align_genes,
    check_blastn,
    mask_sequence,
    write_fasta,
)
from crossgene.models import BedRegion, GeneRecord


def _make_gene(name="GENE_A", chrom="chr1", start=1000, end=2000, strand="+") -> GeneRecord:
    return GeneRecord(
        name=name, gene_id="G001", chrom=chrom,
        start=start, end=end, strand=strand,
        sequence="A" * (end - start),
    )


class TestCheckBlastn:
    def test_both_found(self):
        with patch("crossgene.blastn.shutil.which", return_value="/usr/bin/blastn"):
            check_blastn()

    @pytest.mark.parametrize("missing", ["blastn", "makeblastdb"])
    def test_tool_not_found(self, missing):
        def mock_which(name):
            return None if name == missing else f"/usr/bin/{name}"
        with patch("crossgene.blastn.shutil.which", side_effect=mock_which):
            with pytest.raises(RuntimeError, match=f"{missing} not found"):
                check_blastn()


class TestBuildBlastnCommand:
    @pytest.mark.parametrize("params,expected_word_size,expected_evalue", [
        (BlastParams(sensitivity=1), "7", "10"),
        (BlastParams(sensitivity=2), "7", "1"),
        (BlastParams(sensitivity=3), "11", "0.01"),
        (BlastParams(sensitivity=4), "13", "0.001"),
        (BlastParams(sensitivity=5), "15", "0.001"),
    ])
    def test_preset_params(self, tmp_path, params, expected_word_size, expected_evalue):
        cmd = _build_blastn_command(tmp_path / "q.fa", tmp_path / "db", tmp_path / "out.tsv", params)
        assert cmd[0] == "blastn"
        assert cmd[cmd.index("-word_size") + 1] == expected_word_size
        assert cmd[cmd.index("-evalue") + 1] == expected_evalue

    def test_max_secondary(self, tmp_path):
        params = BlastParams(max_secondary=5)
        cmd = _build_blastn_command(tmp_path / "q.fa", tmp_path / "db", tmp_path / "out.tsv", params)
        assert cmd[cmd.index("-max_target_seqs") + 1] == "5"


class TestMakeBlastDb:
    def test_calls_makeblastdb(self, tmp_path):
        target = tmp_path / "target.fa"
        target.write_text(">gene\nACGT\n")
        mock_result = MagicMock(returncode=0, stderr="")
        with patch("crossgene.blastn.subprocess.run", return_value=mock_result) as mock_run:
            db_path, _ = _make_blast_db(target)
            assert mock_run.call_args[0][0][0] == "makeblastdb"
            assert db_path == target

    def test_raises_on_failure(self, tmp_path):
        target = tmp_path / "target.fa"
        target.write_text(">gene\nACGT\n")
        with patch("crossgene.blastn.subprocess.run", return_value=MagicMock(returncode=1, stderr="err")):
            with pytest.raises(RuntimeError, match="makeblastdb failed"):
                _make_blast_db(target)


class TestWriteFastaAndMask:
    def test_writes_fasta(self):
        gene = _make_gene()
        path = write_fasta(gene, prefix="test")
        content = path.read_text()
        assert content.startswith(">GENE_A\n")
        assert len(content.split("\n")[1]) == 1000
        path.unlink()

    @pytest.mark.parametrize("regions,expected", [
        ([], "ACGTACGT"),
        ([BedRegion(chrom="chr1", start=1002, end=1006)], "ACNNNNGT"),
        ([BedRegion(chrom="chr1", start=1000, end=1002),
          BedRegion(chrom="chr1", start=1006, end=1008)], "NNGTACNN"),
    ])
    def test_mask_sequence(self, regions, expected):
        seq = "ACGTACGT"
        result = mask_sequence(seq, regions, gene_start=1000)
        assert result == expected


class TestAlignGenes:
    def test_full_pipeline(self):
        query = _make_gene("GENE_A")
        target = _make_gene("GENE_B", chrom="chr2", start=5000, end=6000)
        mock_result = MagicMock(returncode=0, stderr="")

        with patch("crossgene.blastn.subprocess.run", return_value=mock_result) as mock_run:
            tsv_path, temp_files = align_genes(query, target, BlastParams())
            assert mock_run.call_count == 2
            assert mock_run.call_args_list[0][0][0][0] == "makeblastdb"
            assert mock_run.call_args_list[1][0][0][0] == "blastn"
            for f in temp_files:
                if f.exists():
                    f.unlink()

    def test_with_blacklist(self):
        query = _make_gene("GENE_A")
        target = _make_gene("GENE_B", chrom="chr2", start=5000, end=6000)
        blacklist = [BedRegion(chrom="chr1", start=1100, end=1200)]
        mock_result = MagicMock(returncode=0, stderr="")

        with patch("crossgene.blastn.subprocess.run", return_value=mock_result):
            _, temp_files = align_genes(query, target, BlastParams(), blacklist_regions=blacklist)
            seq = temp_files[0].read_text().strip().split("\n")[1]
            assert seq[100:200] == "N" * 100
            assert seq[0] == "A"
            for f in temp_files:
                if f.exists():
                    f.unlink()
