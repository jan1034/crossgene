#!/bin/bash
set -euo pipefail

REFS="references"
GENOME="${REFS}/mus_musculus.101.mainChr.fa"
GTF="${REFS}/mus_musculus.101.genes.gtf"
ANN_GTF="${REFS}/mus_musculus.101.mainChr.gtf"
CHROM="${REFS}/mus_musculus.101.chrom.sizes"
SINES="${REFS}/mm10_sines.bed"
FLANKING=5000

COMMON="--genome ${GENOME} --genes-gtf ${GTF} --annotation-gtf ${ANN_GTF} --chrom-sizes ${CHROM} --flanking ${FLANKING}"
GENES="--gene-a Actg1 --gene-b Actg2"

# --- All sensitivity x stringency combinations (no blacklist) ---
for sens in 1 2 3 4 5; do
    for str in 1 2 3 4 5; do
        echo "=== sensitivity=${sens} stringency=${str} (no blacklist) ==="
        crossgene ${GENES} ${COMMON} --bed ${SINES} \
            --sensitivity ${sens} --stringency ${str} \
            --outdir output/actg1-actg2/sens${sens}_str${str}/
    done
done

# --- All sensitivity x stringency combinations (SINEs blacklisted) ---
for sens in 1 2 3 4 5; do
    for str in 1 2 3 4 5; do
        echo "=== sensitivity=${sens} stringency=${str} (SINEs blacklisted) ==="
        crossgene ${GENES} ${COMMON} --blacklist ${SINES} --bed ${SINES} \
            --sensitivity ${sens} --stringency ${str} \
            --outdir output/actg1-actg2-no-sines/sens${sens}_str${str}/
    done
done
