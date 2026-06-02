"""BED9 export for alignment hits with numbered cross-references."""

from __future__ import annotations

import logging
import os
from typing import Callable

from crossgene.models import AlignmentHit, GeneRecord

logger = logging.getLogger("crossgene")

COLOR_SENSE = "70,130,180"       # steelblue — same-sense (+)
COLOR_ANTISENSE = "178,34,34"    # firebrick — antisense (-)


def _write_bed_file(
    path: str, track_name: str, track_desc: str,
    sorted_hits: list[AlignmentHit],
    get_coords: Callable[[AlignmentHit], tuple[str, int, int]],
    primary_prefix: str, secondary_prefix: str,
    primary_num: dict[int, int], secondary_num: dict[int, int],
) -> None:
    with open(path, "w") as f:
        f.write(f'track name="{track_name}" description="{track_desc}" itemRgb=on\n')
        for h in sorted_hits:
            pn = primary_num[id(h)]
            sn = secondary_num[id(h)]
            name = f"{primary_prefix}{pn}-{secondary_prefix}{sn}"
            score = min(int(h.identity * 1000), 1000)
            rgb = COLOR_SENSE if h.strand == "+" else COLOR_ANTISENSE
            chrom, start, end = get_coords(h)
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{h.strand}\t{start}\t{end}\t{rgb}\n")


def write_hit_beds(
    hits: list[AlignmentHit],
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    direction_label: str,       # "AtoB" or "BtoA"
    output_dir: str,
    primary_only: bool = False,
) -> tuple[str, str]:
    """Write two BED9 files (query-side, target-side) with numbered cross-refs.

    Returns paths to (query_bed, target_bed).
    """
    # Filter
    if primary_only:
        hits = [h for h in hits if h.is_primary]

    # Determine prefixes
    if direction_label == "AtoB":
        q_prefix, t_prefix = "A", "B"
    else:
        q_prefix, t_prefix = "B", "A"

    # Assign query-side numbers (sort by query_start asc, tiebreak identity desc)
    query_sorted = sorted(hits, key=lambda h: (h.query_start, -h.identity))
    query_num = {id(h): i + 1 for i, h in enumerate(query_sorted)}

    # Assign target-side numbers (sort by target_start asc, tiebreak identity desc)
    target_sorted = sorted(hits, key=lambda h: (h.target_start, -h.identity))
    target_num = {id(h): i + 1 for i, h in enumerate(target_sorted)}

    # File paths
    q_bed_path = os.path.join(output_dir, f"{query_gene.name}_hits_{direction_label}.bed")
    t_bed_path = os.path.join(output_dir, f"{target_gene.name}_hits_{direction_label}.bed")

    _write_bed_file(
        q_bed_path,
        f"{query_gene.name}_hits_{direction_label}",
        f"{query_gene.name} hits → {target_gene.name}",
        query_sorted,
        get_coords=lambda h: (h.query_chrom, h.query_start, h.query_end),
        primary_prefix=q_prefix, secondary_prefix=t_prefix,
        primary_num=query_num, secondary_num=target_num,
    )

    _write_bed_file(
        t_bed_path,
        f"{target_gene.name}_hits_{direction_label}",
        f"{query_gene.name} hits → {target_gene.name}",
        target_sorted,
        get_coords=lambda h: (h.target_chrom, h.target_start, h.target_end),
        primary_prefix=t_prefix, secondary_prefix=q_prefix,
        primary_num=target_num, secondary_num=query_num,
    )

    logger.info("Wrote hit BED: %s (%d hits)", q_bed_path, len(hits))
    logger.info("Wrote hit BED: %s (%d hits)", t_bed_path, len(hits))

    return q_bed_path, t_bed_path
