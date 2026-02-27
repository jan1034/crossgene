"""CLI entry point for compare-genes."""

import click


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
    click.echo(f"compare-genes v0.1.0 — not yet implemented")
    click.echo(f"Gene A: {gene_a}, Gene B: {gene_b}")
