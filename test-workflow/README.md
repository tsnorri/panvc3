# Test workflow for PanVC 3

This directory contains a simple test workflow that uses PanVC 3 with [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/). The provided input consists of simulated data. The workflow can be run with Snakemake with e.g. `./run-snakemake.sh`.

## Requirements

* [Snakemake](https://snakemake.github.io/) 7.22.0 or newer.

## Inputs

* *genome/genome.fa.gz* a (very) simple test genome generated with http://www.faculty.ucr.edu/~mmaduro/random.htm.
* *merged.vcf.gz* variants generated with [Mutation-Simulator](https://pypi.org/project/Mutation-Simulator) and merged with [bcftools](https://samtools.github.io/bcftools/).
* *index-input/sequences* contains predicted haplotype sequences generated with [vcf2multialign](https://github.com/tsnorri/vcf2multialign) from *genome/genome.fa.gz* and *merged.vcf.gz*.
* *reads* contains simulated reads produced with [Mason](https://www.seqan.de/apps/mason.html) from *genome/genome.fa.gz* and *merged.vcf.gz*.

## Outputs

The aligned reads projected to the test genome co-ordinates will be written to `alignments/alignments.mapq-recalculated.sam.gz`.

## Notes

Mason was run with `mason_simulator -ir genome.fa -n 1000 -iv merged.vcf -o left.fq -or right.fq --verbose --seed 21 --num-threads 10 -oa alignment.bam`
