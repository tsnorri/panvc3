#!/bin/bash

set -euxo pipefail

if [ -f local.mk ]
then
	echo "ERROR: local.mk already exists."
	exit 1
fi

echo "Generating local.mk"
m4 -D CONDA_PREFIX="${PREFIX}" conda/local.mk.m4 > local.mk

echo "Running make"
make -j ${CPU_COUNT} lib/libkqueue/build/libkqueue.a
make -j ${CPU_COUNT} all

echo "Copying build products"
dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp alignment-statistics/alignment_statistics					"${dst_bin}/panvc3_alignments_statistics"
cp convert-bed-positions/convert_bed_positions					"${dst_bin}/panvc3_convert_bed_positions"
cp count-supporting-reads/count_supporting_reads				"${dst_bin}/panvc3_count_supporting_reads"
cp count-supporting-reads/calculate_reference_bias.py			"${dst_bin}/panvc3_calculate_reference_bias.py"
cp index-msa/index_msa											"${dst_bin}/panvc3_index_msa"
cp project-alignments/project_alignments						"${dst_bin}/panvc3_project_alignments"
cp process-alignments/process_alignments						"${dst_bin}/panvc3_process_alignments"
cp recalculate-mapq/recalculate_mapq							"${dst_bin}/panvc3_recalculate_mapq"
cp rewrite-cigar/rewrite_cigar									"${dst_bin}/panvc3_rewrite_cigar"
cp split-alignments-by-reference/split_alignments_by_reference	"${dst_bin}/panvc3_split_alignments_by_reference"
cp subset-alignments/subset_alignments							"${dst_bin}/panvc3_subset_alignments"
