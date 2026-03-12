"""BLASTN alignment wrapper for gene-vs-gene comparison."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from crossgene.models import BedRegion, GeneRecord
from crossgene.presets import get_sensitivity_preset

logger = logging.getLogger(__name__)


@dataclass
class BlastParams:
    """Parameters for BLASTN alignment."""

    max_secondary: int = 10  # maps to max_target_seqs / max_hsps
    sensitivity: int = 3  # 1 (most sensitive) to 5 (strictest)


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


def write_fasta(gene: GeneRecord, prefix: str = "gene") -> Path:
    """Write a gene's sequence to a temporary FASTA file."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", prefix=f"{prefix}_{gene.name}_", delete=False
    )
    tmp.write(f">{gene.name}\n{gene.sequence}\n")
    tmp.close()
    return Path(tmp.name)


def mask_sequence(
    sequence: str, regions: list[BedRegion], gene_start: int,
) -> str:
    """Replace blacklisted regions with N's in the gene sequence.

    Regions are in genomic coordinates, already filtered/clipped to gene bounds.
    The sequence is in genomic orientation, so local position 0 = gene_start.

    Args:
        sequence: Gene sequence string (genomic orientation).
        regions: BedRegions already filtered/clipped to gene bounds.
        gene_start: Genomic start of the gene (0-based).

    Returns:
        Masked sequence with N's replacing blacklisted bases.
    """
    seq = list(sequence)
    for region in regions:
        local_start = region.start - gene_start
        local_end = region.end - gene_start
        for i in range(max(0, local_start), min(len(seq), local_end)):
            seq[i] = 'N'
    return ''.join(seq)


def _build_blastn_command(
    query_fasta: Path, db_path: Path, output_tsv: Path, params: BlastParams
) -> list[str]:
    """Build the blastn command line with tabular output."""
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-out", str(output_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore qlen slen sstrand",
        "-max_target_seqs", str(params.max_secondary),
        "-dust", "yes",
        "-soft_masking", "false",
    ]

    preset = get_sensitivity_preset(params.sensitivity)
    cmd.extend([
        "-word_size", str(preset["word_size"]),
        "-reward", str(preset["reward"]),
        "-penalty", str(preset["penalty"]),
        "-gapopen", str(preset["gapopen"]),
        "-gapextend", str(preset["gapextend"]),
        "-evalue", str(preset["evalue"]),
    ])

    return cmd


def _make_blast_db(target_fasta: Path) -> tuple[Path, list[Path]]:
    """Run makeblastdb on the target FASTA.

    Returns (db_path, list_of_db_files_to_clean).
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

    db_extensions = [".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"]
    db_files = [
        Path(str(target_fasta) + ext)
        for ext in db_extensions
        if Path(str(target_fasta) + ext).exists()
    ]

    return target_fasta, db_files


def align_genes(
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    params: BlastParams,
    blacklist_regions: list[BedRegion] | None = None,
) -> tuple[Path, list[Path]]:
    """Run BLASTN gene-vs-gene alignment.

    Writes query and target to temp FASTAs, builds a BLAST DB from the target,
    and runs BLASTN with the query against it.

    If blacklist_regions are provided, the query sequence is masked with N's
    at those positions before alignment.

    Returns (tsv_path, temp_files_to_clean).
    """
    check_blastn()

    temp_files: list[Path] = []

    # Write query FASTA (with optional blacklist masking)
    if blacklist_regions:
        masked_seq = mask_sequence(
            query_gene.sequence, blacklist_regions, query_gene.start,
        )
        masked_gene = replace(query_gene, sequence=masked_seq)
        query_fasta = write_fasta(masked_gene, prefix="query")
        n_count = sum(1 for a, b in zip(query_gene.sequence, masked_seq) if a != b)
        logger.info("Blacklist: masked %d bases in %s query sequence", n_count, query_gene.name)
    else:
        query_fasta = write_fasta(query_gene, prefix="query")
    temp_files.append(query_fasta)

    # Write target FASTA and build DB
    target_fasta = write_fasta(target_gene, prefix="target")
    temp_files.append(target_fasta)

    logger.debug("Building BLAST database from %s", target_fasta)
    db_path, db_files = _make_blast_db(target_fasta)
    temp_files.extend(db_files)

    # Run BLASTN
    output_tsv = Path(tempfile.mktemp(suffix=".tsv", prefix="blast_"))
    cmd = _build_blastn_command(query_fasta, db_path, output_tsv, params)
    temp_files.append(output_tsv)

    logger.debug("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("blastn: %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"blastn failed (exit {result.returncode}):\n{result.stderr}"
        )

    return output_tsv, temp_files
