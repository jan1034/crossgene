"""CLI entry point for compare-genes."""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import click

from compare_genes import __version__
from compare_genes.align import AlignParams, align_fragments, check_minimap2, write_target_fasta
from compare_genes.bigwig import write_bigwig
from compare_genes.fragment import generate_fragments
from compare_genes.gene_extractor import extract_sequence, lookup_gene
from compare_genes.parser import parse_paf
from compare_genes.scores import compute_scores
from compare_genes.tsv_writer import write_tsv
from compare_genes.visualize import create_circlize_plot

logger = logging.getLogger("compare_genes")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def _parse_formats(output_formats: str) -> set[str]:
    valid = {"bigwig", "tsv", "plot"}
    formats = {f.strip().lower() for f in output_formats.split(",")}
    unknown = formats - valid
    if unknown:
        raise click.BadParameter(
            f"Unknown format(s): {', '.join(unknown)}. Valid: {', '.join(sorted(valid))}"
        )
    return formats


def _run_direction(
    query_gene, target_gene, fragment_size, step_size, min_quality,
    align_params, chrom_sizes, outdir, formats, direction_label,
) -> tuple[list, list[Path]]:
    """Run one direction of the comparison pipeline. Returns (hits, temp_files)."""
    temp_files: list[Path] = []

    logger.info(
        "Direction %s: fragmenting %s (%d bp, step=%d)",
        direction_label, query_gene.name,
        len(query_gene.sequence), step_size,
    )
    frags_path = generate_fragments(query_gene, fragment_size, step_size)
    temp_files.append(frags_path)

    target_path = write_target_fasta(target_gene)
    temp_files.append(target_path)

    logger.info("Aligning fragments to %s...", target_gene.name)
    paf_path = align_fragments(frags_path, target_path, align_params)
    temp_files.append(paf_path)

    hits = parse_paf(
        paf_path, query_gene, target_gene,
        fragment_size, step_size, min_quality, direction_label,
    )
    logger.info("Found %d hits passing quality filter", len(hits))

    if not hits:
        logger.warning("No alignments found for %s", direction_label)

    # Output file naming: QUERY_vs_TARGET
    prefix = f"{query_gene.name}_vs_{target_gene.name}"

    if "bigwig" in formats:
        scores = compute_scores(hits, query_gene, fragment_size, step_size)
        bw_path = os.path.join(outdir, f"{prefix}.mappability.bw")
        write_bigwig(scores, query_gene.chrom, query_gene.start, chrom_sizes, bw_path)
        logger.info("Wrote BigWig: %s", bw_path)

    if "tsv" in formats:
        tsv_path = os.path.join(outdir, f"{prefix}.hits.tsv")
        write_tsv(hits, tsv_path)
        logger.info("Wrote TSV: %s (%d rows)", tsv_path, len(hits))

    return hits, temp_files


@click.command()
@click.option("--gene-a", required=True, help="First gene name (e.g., BRCA1)")
@click.option("--gene-b", required=True, help="Second gene name (e.g., BRCA2)")
@click.option("--fragment-size", default=50, show_default=True, help="Fragment length in bp")
@click.option("--step-size", default=1, show_default=True, help="Step between fragments in bp")
@click.option("--min-quality", default=30, show_default=True, help="Minimum identity (0-100) to report a hit")
@click.option("--max-secondary", default=10, show_default=True, help="Max secondary alignments per fragment")
@click.option("--genome", default="references/homo_sapiens.109.mainChr.fa", show_default=True, help="Path to genome FASTA")
@click.option("--gtf", default="references/homo_sapiens.109.genes.gtf", show_default=True, help="Path to gene annotation GTF")
@click.option("--chrom-sizes", default="references/homo_sapiens.109.chrom.sizes", show_default=True, help="Chromosome sizes file")
@click.option("--outdir", default=".", show_default=True, help="Output directory")
@click.option("--output-formats", default="bigwig,tsv,plot", show_default=True, help="Comma-separated: bigwig,tsv,plot")
@click.option("--minimap2-preset", default="auto", show_default=True, help="minimap2 preset (auto selects based on fragment size)")
@click.option("--sensitive", is_flag=True, default=False, help="Tune minimap2 for short fragments (<=50bp)")
@click.option("--verbose", "-v", is_flag=True, default=False, help="Enable debug logging")
def main(gene_a, gene_b, fragment_size, step_size, min_quality, max_secondary,
         genome, gtf, chrom_sizes, outdir, output_formats, minimap2_preset,
         sensitive, verbose):
    """Compare two gene sequences by fragment alignment.

    Fragments one gene, aligns to another using minimap2, and produces
    similarity profiles (BigWig), hit tables (TSV), and circular visualizations.
    Runs bidirectional comparison (A->B and B->A).
    """
    t0 = time.time()
    _setup_logging(verbose)
    logger.info("compare-genes v%s", __version__)

    # Parse and validate output formats
    formats = _parse_formats(output_formats)

    # Validate inputs
    try:
        check_minimap2()
    except RuntimeError as e:
        raise click.ClickException(str(e))

    for path, label in [(genome, "Genome FASTA"), (gtf, "GTF"), (chrom_sizes, "Chrom sizes")]:
        if not os.path.exists(path):
            raise click.ClickException(f"{label} not found: {path}")

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Extract genes
    logger.info("Looking up genes...")
    try:
        rec_a = lookup_gene(gene_a, gtf)
        rec_b = lookup_gene(gene_b, gtf)
    except ValueError as e:
        raise click.ClickException(str(e))

    logger.info(
        "Gene A: %s %s:%d-%d (%s)",
        rec_a.name, rec_a.chrom, rec_a.start, rec_a.end, rec_a.strand,
    )
    logger.info(
        "Gene B: %s %s:%d-%d (%s)",
        rec_b.name, rec_b.chrom, rec_b.start, rec_b.end, rec_b.strand,
    )

    rec_a = extract_sequence(rec_a, genome)
    rec_b = extract_sequence(rec_b, genome)

    align_params = AlignParams(
        fragment_size=fragment_size,
        max_secondary=max_secondary,
        minimap2_preset=minimap2_preset,
        sensitive=sensitive,
    )

    all_temp_files: list[Path] = []

    # Direction A→B
    hits_ab, temps = _run_direction(
        rec_a, rec_b, fragment_size, step_size, min_quality,
        align_params, chrom_sizes, outdir, formats, "A→B",
    )
    all_temp_files.extend(temps)

    # Direction B→A
    hits_ba, temps = _run_direction(
        rec_b, rec_a, fragment_size, step_size, min_quality,
        align_params, chrom_sizes, outdir, formats, "B→A",
    )
    all_temp_files.extend(temps)

    # Visualization (A→B hits only, per architecture)
    if "plot" in formats:
        plot_path = os.path.join(outdir, f"{rec_a.name}_vs_{rec_b.name}.circlize.pdf")
        logger.info("Creating circlize plot...")
        create_circlize_plot(hits_ab, rec_a, rec_b, plot_path)

    # Clean up temp files
    for tmp in all_temp_files:
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass

    elapsed = time.time() - t0
    logger.info("Done in %.1fs. A→B: %d hits, B→A: %d hits", elapsed, len(hits_ab), len(hits_ba))
