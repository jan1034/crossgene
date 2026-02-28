"""minimap2 alignment wrapper."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from compare_genes.models import GeneRecord

logger = logging.getLogger(__name__)


@dataclass
class AlignParams:
    """Parameters for minimap2 alignment."""

    fragment_size: int = 50
    max_secondary: int = 10
    minimap2_preset: str = "auto"  # "auto", "sr", "map-ont", etc.
    sensitive: bool = False


def check_minimap2() -> None:
    """Verify minimap2 is installed and accessible.

    Raises RuntimeError with install instructions if not found.
    """
    if shutil.which("minimap2") is None:
        raise RuntimeError(
            "minimap2 not found on PATH. Install it with:\n"
            "  conda install -c bioconda minimap2\n"
            "or see https://github.com/lh3/minimap2"
        )


def _select_preset(params: AlignParams) -> list[str]:
    """Select minimap2 preset flags based on parameters."""
    if params.minimap2_preset != "auto":
        return ["-x", params.minimap2_preset]
    if params.fragment_size <= 200:
        return ["-x", "sr"]
    return []


def _build_command(
    fragment_fasta: Path,
    target_fasta: Path,
    output_paf: Path,
    params: AlignParams,
) -> list[str]:
    """Build the minimap2 command line."""
    cmd = [
        "minimap2",
        "-c",           # output CIGAR in PAF
        "--eqx",        # use =/X in CIGAR instead of M
        "--secondary=yes",
        "-N", str(params.max_secondary),
    ]

    cmd.extend(_select_preset(params))

    if params.sensitive:
        cmd.extend(["-k", "11", "-w", "5"])

    cmd.extend([
        "-o", str(output_paf),
        str(target_fasta),
        str(fragment_fasta),
    ])

    return cmd


def write_target_fasta(gene: GeneRecord) -> Path:
    """Write a gene's sequence to a temporary FASTA file for use as minimap2 reference."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", prefix=f"target_{gene.name}_", delete=False
    )
    tmp.write(f">{gene.name}\n{gene.sequence}\n")
    tmp.close()
    return Path(tmp.name)


def align_fragments(
    fragment_fasta: Path, target_fasta: Path, params: AlignParams
) -> Path:
    """Align fragment FASTA against a target FASTA using minimap2.

    Returns the path to the output PAF file.
    Raises RuntimeError if minimap2 fails.
    """
    check_minimap2()

    output_paf = Path(tempfile.mktemp(suffix=".paf", prefix="align_"))
    cmd = _build_command(fragment_fasta, target_fasta, output_paf, params)

    logger.debug("Running: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("minimap2: %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed (exit {result.returncode}):\n{result.stderr}"
        )

    return output_paf
