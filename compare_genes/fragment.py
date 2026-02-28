"""Rolling-window fragmentation of gene sequences."""

from __future__ import annotations

import tempfile
from pathlib import Path

from compare_genes.models import GeneRecord


def generate_fragments(
    gene: GeneRecord, fragment_size: int, step_size: int
) -> Path:
    """Fragment a gene's sense sequence using a rolling window.

    Writes fragments to a temporary FASTA file with headers ``>gene_name:idx``
    where idx is a 0-based sequential index.

    Edge cases:
    - Sequence shorter than fragment_size: single fragment = full sequence.
    - Last fragment shorter than fragment_size: included if its length
      >= fragment_size / 2, otherwise skipped.

    Returns the path to the temporary FASTA file.
    """
    seq = gene.sequence
    seq_len = len(seq)

    if seq_len == 0:
        raise ValueError(f"Gene '{gene.name}' has no sequence; call extract_sequence first")

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", prefix=f"frags_{gene.name}_", delete=False
    )

    idx = 0
    if seq_len <= fragment_size:
        # Sequence shorter than or equal to fragment size: single fragment
        tmp.write(f">{gene.name}:{idx}\n{seq}\n")
    else:
        pos = 0
        while pos + fragment_size <= seq_len:
            frag = seq[pos : pos + fragment_size]
            tmp.write(f">{gene.name}:{idx}\n{frag}\n")
            idx += 1
            pos += step_size

        # Handle trailing partial fragment
        if pos < seq_len:
            remainder = seq[pos:]
            if len(remainder) >= fragment_size / 2:
                tmp.write(f">{gene.name}:{idx}\n{remainder}\n")

    tmp.close()
    return Path(tmp.name)


def fragment_index_to_genomic(
    gene: GeneRecord, idx: int, fragment_size: int, step_size: int
) -> tuple[int, int]:
    """Convert a fragment index to genomic coordinates (0-based half-open).

    For plus-strand genes:
        sense_start = idx * step_size
        genomic = [gene.start + sense_start, gene.start + sense_start + frag_len)

    For minus-strand genes (sense sequence is reverse-complemented):
        sense_start = idx * step_size
        genomic = [gene.end - sense_start - frag_len, gene.end - sense_start)

    The returned fragment length is capped at the actual sequence length
    to handle the trailing partial fragment correctly.
    """
    seq_len = len(gene.sequence)
    sense_start = idx * step_size
    frag_len = min(fragment_size, seq_len - sense_start)

    if gene.strand == "-":
        genomic_end = gene.end - sense_start
        genomic_start = genomic_end - frag_len
    else:
        genomic_start = gene.start + sense_start
        genomic_end = genomic_start + frag_len

    return genomic_start, genomic_end
