# Copyright (c) 2022 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import sys


class Record(object):
	def __init__(self, pos, var_id, ref, alt):
		self.pos = pos
		self.var_id = var_id
		self.ref = ref
		self.alt = alt
		self.support = []

	def var_as_tuple(self):
		return (self.var_id, self.ref, self.alt, self.support)
	
	def support_string(self):
		return " ".join(map(lambda x: f"{x[0]}:{x[1]}", self.support))


class RecordGenerator(object):
	def __init__(self, fp):
		self.fp = fp
		self.is_first = True
		self.current_record = None
	
	def parse(self):
		for line_ in self.fp:
			line = line_.rstrip("\n")
			fields = line.split("\t")
			rec_type = fields[0]

			if 'V' == rec_type:
				if not self.is_first:
					self.current_record.support.sort()
					yield self.current_record
				else:
					self.is_first = False

				pos = int(fields[2])
				var_id, ref, alt = fields[3:6]
				self.current_record = Record(pos, var_id, ref, alt)
			elif 'R' == rec_type:
				count = int(fields[1])
				seq = fields[2]
				self.current_record.support.append((seq, count))

		if not self.is_first:
			self.current_record.support.sort()
			yield self.current_record


def cmp(lhs, rhs):
	return (lhs > rhs) - (lhs < rhs)


def list_merge(lhsc, rhsc, key_fn):
	lhs = iter(lhsc)
	rhs = iter(rhsc)

	ll = next(lhs, None)
	rr = next(rhs, None)

	while not ((ll is None) or (rr is None)):
		res = cmp(key_fn(ll), key_fn(rr))

		if -1 == res:
			yield (ll, None)
			ll = next(lhs, None)
		elif 1 == res:
			yield (None, rr)
			rr = next(rhs, None)
		else:
			yield (ll, rr)
			ll = next(lhs, None)
			rr = next(rhs, None)
	
	if ll is None and rr is None:
		return

	if ll is None:
		yield (None, rr)
		yield from map(lambda x: (None, x), rhs)
		return

	if rr is None:
		yield (ll, None)
		yield from map(lambda x: (x, None), lhs)
		return


def coalesce(*values):
	return next(filter(lambda x: x is not None, values), None)


class Diff(object):

	def parse_one(self, fp):
		gen = RecordGenerator(fp)
		res = list(gen.parse())
		res.sort(key = lambda x: x.var_id)
		return res

	def parse(self, lhs_fp, rhs_fp):
		lhs_recs = self.parse_one(lhs_fp)
		rhs_recs = self.parse_one(rhs_fp)
		paired = list_merge(lhs_recs, rhs_recs, lambda x: x.var_id)
		diff = list(filter(lambda x: x[0] is None or x[1] is None or x[0].var_as_tuple() != x[1].var_as_tuple(), paired))
		diff.sort(key = lambda x: coalesce(*list(x)).pos)

		for lhs, rhs in diff:
			if lhs is None:
				sys.stdout.write(f"{rhs.var_id}\t\t{rhs.ref}\t\t{rhs.alt}\t\t{rhs.support_string()}\n")
				continue

			if rhs is None:
				sys.stdout.write(f"{lhs.var_id}\t{lhs.ref}\t\t{lhs.alt}\t\t{lhs.support_string()}\t\n")

			assert(lhs.var_id == rhs.var_id)
			sys.stdout.write(f"{lhs.var_id}\t{lhs.ref}\t{rhs.ref}\t{lhs.alt}\t{rhs.alt}\t{lhs.support_string()}\t{rhs.support_string()}\n")

		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'Compare two sets of supporting reads by read IDs')
	parser.add_argument('lhs', type = argparse.FileType('r'))
	parser.add_argument('rhs', type = argparse.FileType('r'))
	args = parser.parse_args()

	diff = Diff()
	diff.parse(args.lhs, args.rhs)
