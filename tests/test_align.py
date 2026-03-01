"""Tests for align module."""

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from crossgene.align import (
    AlignParams,
    _build_command,
    _select_preset,
    align_fragments,
    check_minimap2,
    write_target_fasta,
)
from crossgene.models import GeneRecord


class TestSelectPreset:
    def test_auto_short_fragments(self):
        params = AlignParams(fragment_size=50)
        assert _select_preset(params) == ["-x", "sr"]

    def test_auto_200bp(self):
        params = AlignParams(fragment_size=200)
        assert _select_preset(params) == ["-x", "sr"]

    def test_auto_long_fragments(self):
        params = AlignParams(fragment_size=500)
        assert _select_preset(params) == []

    def test_user_override(self):
        params = AlignParams(fragment_size=50, minimap2_preset="map-ont")
        assert _select_preset(params) == ["-x", "map-ont"]


class TestBuildCommand:
    def test_basic_command(self, tmp_path):
        frags = tmp_path / "frags.fa"
        target = tmp_path / "target.fa"
        output = tmp_path / "out.paf"
        params = AlignParams(fragment_size=50, max_secondary=10)

        cmd = _build_command(frags, target, output, params)

        assert cmd[0] == "minimap2"
        assert "-c" in cmd
        assert "--eqx" in cmd
        assert "--secondary=yes" in cmd
        assert "-N" in cmd
        assert "10" in cmd
        assert "-x" in cmd
        assert "sr" in cmd
        assert str(target) in cmd
        assert str(frags) in cmd
        assert "-o" in cmd

    def test_sensitive_mode(self, tmp_path):
        frags = tmp_path / "frags.fa"
        target = tmp_path / "target.fa"
        output = tmp_path / "out.paf"
        params = AlignParams(fragment_size=50, sensitive=True)

        cmd = _build_command(frags, target, output, params)

        assert "-k" in cmd
        assert "9" in cmd
        assert "-w" in cmd
        assert "3" in cmd
        assert "-A" in cmd
        assert "-B" in cmd
        assert "-m" in cmd
        assert "-s" in cmd

    def test_divergent_mode(self, tmp_path):
        frags = tmp_path / "frags.fa"
        target = tmp_path / "target.fa"
        output = tmp_path / "out.paf"
        params = AlignParams(fragment_size=50, divergent=True)

        cmd = _build_command(frags, target, output, params)

        assert "-k" in cmd
        assert "7" in cmd
        assert "-w" in cmd
        assert "2" in cmd
        assert "-A" in cmd
        assert "-B" in cmd
        assert "-m" in cmd
        assert "-s" in cmd

    def test_divergent_overrides_sensitive(self, tmp_path):
        """When divergent is set, its params should be used (not sensitive)."""
        frags = tmp_path / "frags.fa"
        target = tmp_path / "target.fa"
        output = tmp_path / "out.paf"
        params = AlignParams(fragment_size=50, divergent=True, sensitive=False)

        cmd = _build_command(frags, target, output, params)
        k_idx = cmd.index("-k")
        assert cmd[k_idx + 1] == "7"

    def test_no_sensitive_by_default(self, tmp_path):
        frags = tmp_path / "frags.fa"
        target = tmp_path / "target.fa"
        output = tmp_path / "out.paf"
        params = AlignParams(fragment_size=50)

        cmd = _build_command(frags, target, output, params)

        assert "-k" not in cmd
        assert "-w" not in cmd


class TestCheckMinimap2:
    def test_found(self):
        # Should not raise if minimap2 is installed
        check_minimap2()

    def test_not_found(self):
        with patch("crossgene.align.shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="minimap2 not found"):
                check_minimap2()


class TestWriteTargetFasta:
    def test_writes_fasta(self):
        gene = GeneRecord(
            name="TEST", gene_id="G001", chrom="chr1",
            start=0, end=100, strand="+", sequence="ACGT" * 25,
        )
        path = write_target_fasta(gene)
        content = path.read_text()
        assert content.startswith(">TEST\n")
        assert "ACGT" * 25 in content
        path.unlink()


class TestAlignFragments:
    def test_mock_subprocess(self, tmp_path):
        """Verify subprocess is called with correct arguments."""
        frags = tmp_path / "frags.fa"
        frags.write_text(">frag:0\nACGT\n")
        target = tmp_path / "target.fa"
        target.write_text(">target\nACGTACGT\n")
        params = AlignParams(fragment_size=50, max_secondary=5)

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""

        with patch("crossgene.align.subprocess.run", return_value=mock_result) as mock_run:
            paf = align_fragments(frags, target, params)

            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "minimap2"
            assert "-N" in call_args
            idx = call_args.index("-N")
            assert call_args[idx + 1] == "5"

            # Clean up
            if paf.exists():
                paf.unlink()

    def test_raises_on_failure(self, tmp_path):
        frags = tmp_path / "frags.fa"
        frags.write_text(">frag:0\nACGT\n")
        target = tmp_path / "target.fa"
        target.write_text(">target\nACGTACGT\n")
        params = AlignParams()

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "some error"

        with patch("crossgene.align.subprocess.run", return_value=mock_result):
            with pytest.raises(RuntimeError, match="minimap2 failed"):
                align_fragments(frags, target, params)
