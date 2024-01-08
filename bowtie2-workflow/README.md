# Workflow for PanVC 3 and Bowtie 2

This directory contains a configurable workflow that uses PanVC 3 with [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/). The example configuration file uses the inputs of the [test workflow](../test-workflow/). Please note that since the graph produced from the example input does not have any possible cut positions, vcf2multialign will use some of the input sequences as the founder sequences.

## Requirements

* Linux on x86-64 with [glibc](https://www.gnu.org/software/libc/) 2.28 or newer to use the binaries available in [Anaconda](https://anaconda.org). (`ldd --version` may be used to check the version of glibc installed with your operating system.)
* [Snakemake](https://snakemake.github.io/) 7.34.0 or newer.

## Configuration

The Snakefile reads the following values from the input configuration.

### Input options

The known variants are read from VCF files in a set of paths. Each path is produced by concatenating `known_variants_prefix`, a chromosome identifier, and `known_variants_suffix`. The chromosome identifiers are listed in `chromosomes`.

| Key                     | Explanation                                                                                                      |
|-------------------------|------------------------------------------------------------------------------------------------------------------|
| `reference_fasta`       | Path to the uncompressed reference FASTA                                                                         |
| `reads_1`               | Input reads as gzipped FASTQ, mate 1                                                                             |
| `reads_2`               | Input reads as gzipped FASTQ, mate 2                                                                             |
| `known_variants_prefix` | Prefix of the path to the uncompressed VCF that contains the known variants for generating the founder sequences |
| `known_variants_suffix` | Suffix of the aforementioned path                                                                                |
| `chromosomes`           | List of chromosomes for known variants                                                                           |

### Founder sequence options

| Key                                  | Explanation                                                           |
|--------------------------------------|-----------------------------------------------------------------------|
| `founder_count`                      | Founder sequence count (in addition to the reference sequence)        |
| `founder_minimum_component_distance` | Minimum component distance between cut positions in the variant graph |

### Output options

| Key            | Explanation                                                    |
|----------------|----------------------------------------------------------------|
| `output_root`  | Output path prefix, may be set to the empty string             |
| `alignment_id` | User-specified (sample) identifier, used as part of file names |

## Running with the example configuration

1. `git clone https://github.com/tsnorri/panvc3.git`
2. `cd panvc3/bowtie2-workflow`
3. Test running Snakemake with e.g. `./run-snakemake.sh`
4. Copy or modify the parameters and run Snakemake without `--dry-run`.

## Outputs

The `all` target causes the input reads to be aligned, the alignemnts to be projected to the specified standard reference and the mapping qualities to be recalculated. The output will be placed to *output_root*/alignments/*alignment_id*.panvc3-bowtie2-f*F*-d*D*.mapq-recalculated.sam.gz where *F* is the number of founder sequences and *D* is the distance between graph components as specified in the configuration.
