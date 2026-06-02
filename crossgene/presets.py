"""Preset configurations for sensitivity and stringency levels."""

from __future__ import annotations

import json
from importlib import resources


def _load_presets() -> dict:
    """Load presets from the bundled JSON file."""
    ref = resources.files("crossgene").joinpath("presets.json")
    return json.loads(ref.read_text(encoding="utf-8"))


def get_sensitivity_preset(level: int) -> dict:
    """Return BLAST parameters for the given sensitivity level (1-5)."""
    if level < 1 or level > 5:
        raise ValueError(f"Sensitivity level must be 1-5, got {level}")
    return _load_presets()["sensitivity"][str(level)]


def get_stringency_preset(level: int) -> dict:
    """Return filter parameters for the given stringency level (1-5).

    Converts string "inf" values to float("inf").
    """
    if level < 1 or level > 5:
        raise ValueError(f"Stringency level must be 1-5, got {level}")
    preset = _load_presets()["stringency"][str(level)]
    for key, value in preset.items():
        if value == "inf":
            preset[key] = float("inf")
    return preset
