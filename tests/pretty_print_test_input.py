# Copyright (c) 2022 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

# Input: test case input for testing alignment_projector or rewrite_cigar().
# Output: pretty-printed test input.

import sys


def output_indices(seq):
	idx = 0
	for cc in seq:
		if '-' == cc:
			sys.stdout.write(' ')
		else:
			sys.stdout.write(str(idx % 10))
			idx += 1


def output_seq(seq):
	sys.stdout.write(seq)


for line in sys.stdin:
	line = line.rstrip("\n")
	if line.startswith("S\t"):
		_, entry_name, lhs_seq, rhs_seq = line.split("\t")

		sys.stdout.write(80 * '=')
		sys.stdout.write('\n')

		sys.stdout.write(f"{entry_name}\n\n")

		output_indices(lhs_seq)
		sys.stdout.write('\n')
		output_seq(lhs_seq)
		sys.stdout.write('\n')

		output_seq(rhs_seq)
		sys.stdout.write('\n')
		output_indices(rhs_seq)
		sys.stdout.write('\n\n')
	elif line.startswith("Q\t"):
		_, query, src_pos, dst_pos, input_cigar, output_cigar, section_name, segment_name = line.split("\t")

		src_pos = int(src_pos)
		dst_pos = int(dst_pos)

		sys.stdout.write(src_pos * ' ')
		sys.stdout.write(query)
		sys.stdout.write('\n')

