# Copyright (c) 2022 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import sys


def cmp(lhs, rhs):
	return (lhs > rhs) - (lhs < rhs)


def merge(lhsc, rhsc, key_fn, merge_fn):
	# Keys are expected to be unique in each of lhsc, rhsc.
	lhs = iter(lhsc)
	rhs = iter(rhsc)

	ll = next(lhs)
	rr = next(rhs)

	while (not (ll is None or rr is None)):
		res = cmp(key_fn(ll), key_fn(rr))
		if -1 == res:
			yield ll
			ll = next(lhs, None)
		elif 1 == res:
			yield rr
			rr = next(rhs, None)
		else:
			ll = merge_fn(ll, rr)
			rr = next(rhs, None)
	
	if ll is None and rr is None:
		return
	
	if ll is None:
		yield rr
		yield from rhs
		return
	
	if rr is None:
		yield ll
		yield from lhs
		return


def record_generator(fp):
	for line_ in fp:
		line = line_.rstrip("\n")
		fields = line.split("\t")
		pos, cov = fields
		pos = int(pos)
		cov = int(cov)
		yield pos, cov


def coverage_sum(lhs_, rhs_):

	# Skip the header.
	next(lhs_)
	next(rhs_)

	lhs = record_generator(lhs_)
	rhs = record_generator(rhs_)

	sys.stdout.write("POSITION\tCOVERAGE\n")
	rng = merge(lhs, rhs, lambda x: x[0], lambda x, y: (x[0], x[1] + y[1]))
	for pos, count in rng:
		sys.stdout.write(f"{pos}\t{count}\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'Merge two sets of coverages')
	parser.add_argument('lhs', type = argparse.FileType('r'))
	parser.add_argument('rhs', type = argparse.FileType('r'))
	args = parser.parse_args()

	coverage_sum(args.lhs, args.rhs)
