"""BLASTN alignment wrapper."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class BlastParams:
    """Parameters for BLASTN alignment."""

    fragment_size: int = 50
    max_secondary: int = 10  # maps to max_target_seqs
    sensitive: bool = False
    divergent: bool = False


def check_blastn() -> None:
    """Verify blastn and makeblastdb are installed and accessible.

    Raises RuntimeError with install instructions if not found.
    """
    for tool in ("blastn", "makeblastdb"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"{tool} not found on PATH. Install BLAST+ with:\n"
                "  conda install -c bioconda blast\n"
                "or see https://blast.ncbi.nlm.nih.gov/doc/blast-help/"
            )


def _build_blastn_command(
    fragment_fasta: Path, db_path: Path, output_tsv: Path, params: BlastParams
) -> list[str]:
    """Build the blastn command line with tabular output."""
    cmd = [
        "blastn",
        "-query", str(fragment_fasta),
        "-db", str(db_path),
        "-out", str(output_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore qlen slen sstrand",
        "-max_target_seqs", str(params.max_secondary),
        "-dust", "no",
        "-soft_masking", "false",
    ]

    if params.divergent:
        cmd.extend([
            "-word_size", "7",
            "-reward", "1",
            "-penalty", "-1",
            "-gapopen", "2",
            "-gapextend", "1",
            "-evalue", "10",
        ])
    elif params.sensitive:
        cmd.extend([
            "-word_size", "7",
            "-reward", "1",
            "-penalty", "-1",
            "-gapopen", "2",
            "-gapextend", "1",
            "-evalue", "1",
        ])
    else:
        cmd.extend([
            "-word_size", "11",
            "-reward", "2",
            "-penalty", "-3",
            "-gapopen", "5",
            "-gapextend", "2",
            "-evalue", "0.1",
        ])

    return cmd


def _make_blast_db(target_fasta: Path) -> tuple[Path, list[Path]]:
    """Run makeblastdb on the target FASTA.

    Returns (db_path, list_of_db_files_to_clean).
    The db_path is the same as target_fasta (BLAST uses it as the base name).
    """
    cmd = [
        "makeblastdb",
        "-in", str(target_fasta),
        "-dbtype", "nucl",
    ]

    logger.debug("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("makeblastdb: %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"makeblastdb failed (exit {result.returncode}):\n{result.stderr}"
        )

    # makeblastdb creates files with extensions .ndb, .nhr, .nin, .not, .nsq, .ntf, .nto
    db_extensions = [".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"]
    db_files = [
        Path(str(target_fasta) + ext)
        for ext in db_extensions
        if Path(str(target_fasta) + ext).exists()
    ]

    return target_fasta, db_files


def align_fragments_blastn(
    fragment_fasta: Path, target_fasta: Path, params: BlastParams
) -> tuple[Path, list[Path]]:
    """Full BLASTN pipeline: makeblastdb -> blastn -> return (tsv_path, temp_files).

    Returns the path to the output TSV and a list of temporary DB files to clean up.
    Raises RuntimeError if any step fails.
    """
    check_blastn()

    # Build BLAST database
    logger.debug("Building BLAST database from %s", target_fasta)
    db_path, db_files = _make_blast_db(target_fasta)

    # Run BLASTN
    output_tsv = Path(tempfile.mktemp(suffix=".tsv", prefix="blast_"))
    cmd = _build_blastn_command(fragment_fasta, db_path, output_tsv, params)

    logger.debug("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("blastn: %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"blastn failed (exit {result.returncode}):\n{result.stderr}"
        )

    return output_tsv, db_files
