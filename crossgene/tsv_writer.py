"""TSV writer for alignment hits."""

from __future__ import annotations

import csv

from crossgene.models import AlignmentHit

COLUMNS = [
    "query_gene", "query_chrom", "frag_start", "frag_end",
    "target_gene", "target_chrom", "hit_start", "hit_end",
    "strand", "identity", "query_coverage", "mapq", "alignment_score",
    "is_primary", "cigar", "direction",
]


def write_tsv(hits: list[AlignmentHit], output_path: str) -> None:
    """Write alignment hits as a TSV file.

    Sorted by frag_start ascending, then alignment_score descending.
    """
    sorted_hits = sorted(hits, key=lambda h: (h.query_start, -h.alignment_score))

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(COLUMNS)

        for h in sorted_hits:
            writer.writerow([
                h.query_gene,
                h.query_chrom,
                h.query_start,
                h.query_end,
                h.target_gene,
                h.target_chrom,
                h.target_start,
                h.target_end,
                h.strand,
                f"{h.identity:.4f}",
                f"{h.query_coverage:.4f}",
                h.mapq,
                h.alignment_score,
                str(h.is_primary).lower(),
                h.cigar,
                h.direction,
            ])
