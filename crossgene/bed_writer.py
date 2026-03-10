"""BED9 export for alignment hits with numbered cross-references."""

from __future__ import annotations

import logging
import os

from crossgene.models import AlignmentHit, GeneRecord

logger = logging.getLogger("crossgene")

COLOR_SENSE = "70,130,180"       # steelblue — same-sense (+)
COLOR_ANTISENSE = "178,34,34"    # firebrick — antisense (-)


def write_hit_beds(
    hits: list[AlignmentHit],
    query_gene: GeneRecord,
    target_gene: GeneRecord,
    direction_label: str,       # "AtoB" or "BtoA"
    output_dir: str,
    primary_only: bool = True,
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

    # Write query BED
    with open(q_bed_path, "w") as f:
        f.write(
            f'track name="{query_gene.name}_hits_{direction_label}" '
            f'description="{query_gene.name} hits → {target_gene.name}" itemRgb=on\n'
        )
        for h in query_sorted:
            qn = query_num[id(h)]
            tn = target_num[id(h)]
            name = f"{q_prefix}{qn}-{t_prefix}{tn}"
            score = min(int(h.identity * 1000), 1000)
            rgb = COLOR_SENSE if h.strand == "+" else COLOR_ANTISENSE
            f.write(
                f"{h.query_chrom}\t{h.query_start}\t{h.query_end}\t{name}\t"
                f"{score}\t{h.strand}\t{h.query_start}\t{h.query_end}\t{rgb}\n"
            )

    # Write target BED
    with open(t_bed_path, "w") as f:
        f.write(
            f'track name="{target_gene.name}_hits_{direction_label}" '
            f'description="{query_gene.name} hits → {target_gene.name}" itemRgb=on\n'
        )
        for h in target_sorted:
            qn = query_num[id(h)]
            tn = target_num[id(h)]
            name = f"{t_prefix}{tn}-{q_prefix}{qn}"
            score = min(int(h.identity * 1000), 1000)
            rgb = COLOR_SENSE if h.strand == "+" else COLOR_ANTISENSE
            f.write(
                f"{h.target_chrom}\t{h.target_start}\t{h.target_end}\t{name}\t"
                f"{score}\t{h.strand}\t{h.target_start}\t{h.target_end}\t{rgb}\n"
            )

    logger.info("Wrote hit BED: %s (%d hits)", q_bed_path, len(hits))
    logger.info("Wrote hit BED: %s (%d hits)", t_bed_path, len(hits))

    return q_bed_path, t_bed_path
