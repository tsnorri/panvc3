# PanVC 3

PanVC 3 is a set of tools to be used as part of a variant calling workflow that uses short reads as its input. The reads are aligned to an index generated from a multiple sequence alignment. A suitable index may be built from founder sequences.

Running a variant calling workflow that utilises PanVC may consist of e.g. the following phases:

- Generating founder sequences from known variants
- Indexing the founder sequences
- Running the read alignment and variant calling workflow

The founder sequences may be generated with [vcf2multialign](https://github.com/tsnorri/vcf2multialign).

## Academic Use

If you use the software in an academic setting, we kindly ask you to cite [Tackling reference bias in genotyping by using founder sequences with PanVC 3](https://doi.org/10.1093/bioadv/vbae027).

```TeX
@article{Norri2024TacklingReferenceBias,
  author = {Norri, Tuukka and Mäkinen, Veli},
  title = {Tackling reference bias in genotyping by using founder sequences with PanVC 3},
  journal = {Bioinformatics Advances},
  volume = {4},
  number = {1},
  pages = {vbae027},
  year = {2024},
  month = {03},
  issn = {2635-0041},
  doi = {10.1093/bioadv/vbae027},
  url = {https://doi.org/10.1093/bioadv/vbae027},
  eprint = {https://academic.oup.com/bioinformaticsadvances/article-pdf/4/1/vbae027/56912765/vbae027.pdf},
}
```

## Running

A simple example workflow and test data are provided in the [test-workflow](test-workflow) subdirectory. The workflow downloads PanVC 3 as well as other required software automatically from [Anaconda](https://anaconda.org). Please see [README.md](test-workflow/README.md) in the subdirectory.

A more complex workflow that uses Bowtie 2 and loads the settings using Snakemake’s configuration (e.g. a YAML file) is in the [bowtie2-workflow](bowtie2-workflow) subdirectory. Please see [README.md](bowtie2-workflow/README.md) in the subdirectory.

## Contents

- [index_msa](index-msa) builds a co-ordinate transformation data structure from a multiple sequence alignment, as well as the sequences as unaligned FASTA to be used as input for a read aligner.
- [project_alignments](project-alignments) uses the co-ordinate transformation data structure to project alignments in BAM or SAM format to well-known co-ordinates, rewrites the CIGAR strings to match the new reference sequence, and realigns parts of the reads if needed.
- [recalculate_mapq](recalculate-mapq) recalculates the mapping qualities of the alignments given as input, taking into account the projected co-ordinate of each alignment.
- [subset_alignments](subset-alignments) subsets the given alignments by some criteria, e.g. selecting the (paired) alignment with the best mapping quality for each read.
- [count_supporting_reads](count-supporting-reads) counts the number of aligned reads that support some known variants. From the output, reference bias can be calculated with [calculate_reference_bias.py](count-supporting-reads/calculate_reference_bias.py).
- [rewrite_cigar](rewrite-cigar) replaces sequence match operations in CIGAR strings (`=` and `X`) with alignment match operations (`M`) and vice-versa.

Please use the `--help` option with each of the tools for usage. See also the [workflow written for the test data](test-workflow/Snakefile).

## Installing

Binaries for Linux on x86-64 are available on [Anaconda](https://anaconda.org). PanVC 3 may be installed with `conda install -c tsnorri -c conda-forge panvc3=v1.0`. [glibc](https://www.gnu.org/software/libc/) 2.28 or newer is required. (`ldd --version` may be used to check the version installed with your operating system.)

## Building

To clone the repository with submodules, please use `git clone --recursive https://github.com/tsnorri/panvc3.git`.

### With [conda-build](https://docs.conda.io/projects/conda-build/en/stable/index.html)

A conda package can be built with conda-build as follows. The build script has been tested with conda-build 3.25.0. [glibc](https://www.gnu.org/software/libc/) 2.28 or newer is required.

1. `cd conda`
2. `./conda-build.sh`

Conda-build will then report the location of the package from which binaries may be extracted.

### By Hand

The following software and libraries are required to build PanVC 3. The tested versions are also listed.

- [Boost 1.82.0](https://www.boost.org)
- [libbz2 1.0.8](https://sourceware.org/bzip2/)
- [Clang C and C++ compilers, version 16.0.6](https://clang.llvm.org)
- [GCC C and C++ compilers and libstdc++, version 12.3.0](https://gcc.gnu.org)
- [GNU Gengetopt 2.23](https://www.gnu.org/software/gengetopt/gengetopt.html)
- [Ragel 6.10](http://www.colm.net/open-source/ragel/)
- [zlib 1.2.13](https://zlib.net)

The following are needed to build [libdispatch](https://apple.github.io/swift-corelibs-libdispatch/) (provided as a Git submodule) on Linux:

- [Cmake 3.26.4](https://cmake.org)
- [Ninja 1.11.1](https://ninja-build.org)

After installing the prerequisites, please do the following:

1. Create a file called `local.mk` in the root of the cloned repository to specify build variables. One of the files [linux-static.local.mk](linux-static.local.mk) and [conda/local.mk.m4](conda/local.mk.m4) may be used as a starting point.
2. Run Make with e.g. `make -j16`.
