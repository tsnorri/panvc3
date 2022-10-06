# Copyright (c) 2022 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import sys
from collections import defaultdict


def histogram(fp, max_length, should_treat_missing_positions_as_zero_coverage):
	# Skip the header.
	next(fp)

	prev_pos_1 = 0
	max_cov = 0
	retval = defaultdict(lambda: 0)
	
	sys.stderr.write("Processing…\n")
	for lineno, line_ in enumerate(fp, start = 1):
		line = line_.rstrip("\n")
		pos, cov = line.split("\t")
		pos = int(pos)
		cov = int(cov)

		retval[cov] += 1
		if should_treat_missing_positions_as_zero_coverage:
			retval[0] += pos - prev_pos_1

		prev_pos_1 = 1 + pos
		max_cov = max(cov, max_cov)

		if 0 == lineno % 1000000:
			sys.stderr.write(f"Processed {lineno} lines…\n")
	
	if should_treat_missing_positions_as_zero_coverage and 0 < max_length:
		retval[0] += 1 + max_length - prev_pos_1
	
	return retval, max_cov


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'Calculate the histogram of coverages.')
	parser.add_argument('--max-length', type = int, default = 0, help = 'Max. sequence length')
	parser.add_argument('--count-missing', action = argparse.BooleanOptionalAction, help = 'Count zeros for missing positions')
	args = parser.parse_args()

	hh, max_cov = histogram(sys.stdin, args.max_length, args.count_missing)

	sys.stdout.write("COVERAGE\tCOUNT\n")
	for i in range(1 + max_cov):
		sys.stdout.write(f"{i}\t{hh[i]}\n")
