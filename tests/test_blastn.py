"""Tests for blastn module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from crossgene.blastn import (
    BlastParams,
    _build_blastn_command,
    _make_blast_db,
    align_fragments_blastn,
    check_blastn,
)


class TestCheckBlastn:
    def test_both_found(self):
        with patch("crossgene.blastn.shutil.which", return_value="/usr/bin/blastn"):
            check_blastn()  # should not raise

    def test_blastn_not_found(self):
        def mock_which(name):
            return None if name == "blastn" else "/usr/bin/makeblastdb"

        with patch("crossgene.blastn.shutil.which", side_effect=mock_which):
            with pytest.raises(RuntimeError, match="blastn not found"):
                check_blastn()

    def test_makeblastdb_not_found(self):
        def mock_which(name):
            return "/usr/bin/blastn" if name == "blastn" else None

        with patch("crossgene.blastn.shutil.which", side_effect=mock_which):
            with pytest.raises(RuntimeError, match="makeblastdb not found"):
                check_blastn()


class TestBuildBlastnCommand:
    def test_default_params(self, tmp_path):
        frags = tmp_path / "frags.fa"
        db = tmp_path / "target"
        output = tmp_path / "out.tsv"
        params = BlastParams()

        cmd = _build_blastn_command(frags, db, output, params)

        assert cmd[0] == "blastn"
        assert "-query" in cmd
        assert str(frags) in cmd
        assert "-db" in cmd
        assert str(db) in cmd
        assert "-out" in cmd
        assert str(output) in cmd
        assert "-max_target_seqs" in cmd
        assert "10" in cmd
        assert "-dust" in cmd
        assert "yes" in cmd

        # Default params
        idx = cmd.index("-word_size")
        assert cmd[idx + 1] == "15"
        idx = cmd.index("-evalue")
        assert cmd[idx + 1] == "1e-3"
        idx = cmd.index("-reward")
        assert cmd[idx + 1] == "1"
        idx = cmd.index("-penalty")
        assert cmd[idx + 1] == "-4"

    def test_sensitive_params(self, tmp_path):
        frags = tmp_path / "frags.fa"
        db = tmp_path / "target"
        output = tmp_path / "out.tsv"
        params = BlastParams(sensitive=True)

        cmd = _build_blastn_command(frags, db, output, params)

        idx = cmd.index("-word_size")
        assert cmd[idx + 1] == "7"
        idx = cmd.index("-evalue")
        assert cmd[idx + 1] == "1"
        idx = cmd.index("-reward")
        assert cmd[idx + 1] == "1"
        idx = cmd.index("-penalty")
        assert cmd[idx + 1] == "-1"
        idx = cmd.index("-gapopen")
        assert cmd[idx + 1] == "2"
        idx = cmd.index("-gapextend")
        assert cmd[idx + 1] == "1"

    def test_divergent_params(self, tmp_path):
        frags = tmp_path / "frags.fa"
        db = tmp_path / "target"
        output = tmp_path / "out.tsv"
        params = BlastParams(divergent=True)

        cmd = _build_blastn_command(frags, db, output, params)

        idx = cmd.index("-word_size")
        assert cmd[idx + 1] == "7"
        idx = cmd.index("-evalue")
        assert cmd[idx + 1] == "10"
        idx = cmd.index("-reward")
        assert cmd[idx + 1] == "1"
        idx = cmd.index("-penalty")
        assert cmd[idx + 1] == "-1"

    def test_max_secondary(self, tmp_path):
        frags = tmp_path / "frags.fa"
        db = tmp_path / "target"
        output = tmp_path / "out.tsv"
        params = BlastParams(max_secondary=5)

        cmd = _build_blastn_command(frags, db, output, params)
        idx = cmd.index("-max_target_seqs")
        assert cmd[idx + 1] == "5"


class TestMakeBlastDb:
    def test_calls_makeblastdb(self, tmp_path):
        target = tmp_path / "target.fa"
        target.write_text(">gene\nACGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""

        with patch("crossgene.blastn.subprocess.run", return_value=mock_result) as mock_run:
            db_path, db_files = _make_blast_db(target)

            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "makeblastdb"
            assert "-in" in call_args
            assert str(target) in call_args
            assert "-dbtype" in call_args
            assert "nucl" in call_args
            assert db_path == target

    def test_raises_on_failure(self, tmp_path):
        target = tmp_path / "target.fa"
        target.write_text(">gene\nACGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "some error"

        with patch("crossgene.blastn.subprocess.run", return_value=mock_result):
            with pytest.raises(RuntimeError, match="makeblastdb failed"):
                _make_blast_db(target)


class TestAlignFragmentsBlastn:
    def test_full_pipeline(self, tmp_path):
        frags = tmp_path / "frags.fa"
        frags.write_text(">gene:0\nACGT\n")
        target = tmp_path / "target.fa"
        target.write_text(">target\nACGTACGT\n")
        params = BlastParams()

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""

        with patch("crossgene.blastn.subprocess.run", return_value=mock_result) as mock_run:
            tsv_path, db_files = align_fragments_blastn(frags, target, params)

            # Should have been called twice: makeblastdb + blastn
            assert mock_run.call_count == 2
            first_call = mock_run.call_args_list[0][0][0]
            assert first_call[0] == "makeblastdb"
            second_call = mock_run.call_args_list[1][0][0]
            assert second_call[0] == "blastn"

            # Clean up
            if tsv_path.exists():
                tsv_path.unlink()

    def test_raises_on_blastn_failure(self, tmp_path):
        frags = tmp_path / "frags.fa"
        frags.write_text(">gene:0\nACGT\n")
        target = tmp_path / "target.fa"
        target.write_text(">target\nACGTACGT\n")
        params = BlastParams()

        # makeblastdb succeeds, blastn fails
        mock_ok = MagicMock()
        mock_ok.returncode = 0
        mock_ok.stderr = ""
        mock_fail = MagicMock()
        mock_fail.returncode = 1
        mock_fail.stderr = "blast error"

        with patch("crossgene.blastn.subprocess.run", side_effect=[mock_ok, mock_fail]):
            with pytest.raises(RuntimeError, match="blastn failed"):
                align_fragments_blastn(frags, target, params)
