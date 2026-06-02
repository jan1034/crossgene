"""Tests for presets module."""

import pytest

from crossgene.presets import get_sensitivity_preset, get_stringency_preset


class TestSensitivityPresets:
    @pytest.mark.parametrize("level", [1, 2, 3, 4, 5])
    def test_all_levels_load(self, level):
        preset = get_sensitivity_preset(level)
        assert "word_size" in preset
        assert "reward" in preset
        assert "penalty" in preset
        assert "gapopen" in preset
        assert "gapextend" in preset
        assert "evalue" in preset

    def test_level_3_matches_moderate(self):
        preset = get_sensitivity_preset(3)
        assert preset["word_size"] == 11
        assert preset["penalty"] == -2
        assert preset["evalue"] == 0.01

    @pytest.mark.parametrize("level", [0, 6, -1])
    def test_invalid_level_raises(self, level):
        with pytest.raises(ValueError, match="Sensitivity level must be 1-5"):
            get_sensitivity_preset(level)


class TestStringencyPresets:
    @pytest.mark.parametrize("level", [1, 2, 3, 4, 5])
    def test_all_levels_load(self, level):
        preset = get_stringency_preset(level)
        assert "min_quality" in preset
        assert "min_mapq" in preset
        assert "min_length" in preset
        assert "min_bitscore" in preset
        assert "max_evalue" in preset

    def test_level_3_matches_defaults(self):
        preset = get_stringency_preset(3)
        assert preset["min_quality"] == 30
        assert preset["min_mapq"] == 0
        assert preset["min_length"] == 25
        assert preset["min_bitscore"] == 0.0
        assert preset["max_evalue"] == float("inf")

    def test_inf_conversion(self):
        preset = get_stringency_preset(1)
        assert preset["max_evalue"] == float("inf")

    @pytest.mark.parametrize("level", [0, 6, -1])
    def test_invalid_level_raises(self, level):
        with pytest.raises(ValueError, match="Stringency level must be 1-5"):
            get_stringency_preset(level)
