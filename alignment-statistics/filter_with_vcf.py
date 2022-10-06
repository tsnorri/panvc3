# Copyright (c) 2022 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

from bisect import bisect_left
import argparse
import functools
import itertools
import sys
import vcfpy


# Two intervals are equivalent iff. they overlap.
class Interval(object):
	def __init__(self, pos, end):
		assert(pos < end)
		self.pos = pos
		self.end = end # Half-open
	
	def __str__(self):
		return f"[{self.pos}, {self.end})"
	
	def __lt__(self, other):
		return self.end <= other.pos
	
	def __gt__(self, other):
		return other.end <= self.pos
	
	def __le__(self, other):
		return not(other < self)
	
	def __ge__(self, other):
		return not(self < other)

	def __eq__(self, other):
		return not(self < other or other < self)
	
	def __ne__(self, other):
		return not(self == other)


def read_bed(bed_fp, chr_id: str):
	""" Build a list of regions, given a BED file."""
	def intervals():
		for line_ in bed_fp:
			line = line_.rstrip("\n")
			fields = line.split("\t")
			chrom = fields[0]
			if chrom != chr_id:
				continue
			chrom_start = int(fields[1])
			chrom_end = int(fields[2])
			yield Interval(chrom_start, chrom_end)
	retval = list(intervals())
	retval.sort()
	return retval


def read_positions(pos_fp):
	for line_ in pos_fp:
		line = line_.rstrip("\n")
		pos, cov = line.split("\t")
		yield int(pos), int(cov)


def filter_vcf_with_chr(vcf_recs, chr_id):
	"""Filter VCF records by chromosome ID."""
	assert vcf_recs is not None
	assert chr_id is not None

	for rec in vcf_recs:
		# Check the chromosome ID.
		if rec.CHROM != chr_id:
			continue

		yield rec


def filter_vcf_with_regions(vcf_recs, regions):
	"""Filter VCF records by regions."""
	assert vcf_recs is not None
	assert regions is not None

	for rec in vcf_recs:
		pos0 = rec.POS - 1
		rec_interval = Interval(pos0, pos0 + len(rec.REF))

		if len(regions) == bisect_left(regions, rec_interval):
			continue

		yield rec


def vcf_position_coverage(vcf_recs, pos_recs):
	"""Find the position records contained in the given VCF records."""
	assert vcf_recs is not None
	assert pos_recs is not None

	pos_rec = next(pos_recs, None)
	def find_next_pos_rec(ii):
		nonlocal pos_recs, pos_rec
		while pos_rec[0] < ii:
			pos_rec = next(pos_recs, None)
			if pos_rec is None:
				return False
		return True
	
	for vcf_rec in vcf_recs:
		vcf_pos0 = vcf_rec.POS - 1
		vcf_end = vcf_pos0 + len(vcf_rec.REF)

		for ii in range(vcf_pos0, vcf_end):
			if find_next_pos_rec(ii):
				if ii < pos_rec[0]:
					yield ii, 0
				else:
					assert ii == pos_rec[0]
					yield pos_rec
			else:
				yield ii, 0


def compose(fns):
	return functools.reduce(lambda f, g: lambda xs: g(f(xs)), fns)


def filter_positions(pos_fp, vcf_path, bed_fp, chr_id):
	regions = None
	if bed_fp is not None:
		regions = read_bed(bed_fp, chr_id)
	
	vcf_reader = vcfpy.Reader.from_path(vcf_path)
	assert vcf_reader is not None
	
	# Header.
	sys.stdout.write(next(pos_fp))
	
	gen_fn = compose(filter(
		lambda x: x is not None,
		[
			(lambda xs: filter_vcf_with_chr(xs, chr_id)) if chr_id is not None else None,
			(lambda xs: filter_vcf_with_regions(xs, regions)) if regions is not None else None,
			lambda xs: vcf_position_coverage(xs, read_positions(pos_fp))
		]
	))
	
	for rec in gen_fn(vcf_reader):
		sys.stdout.write(f"{rec[0]}\t{rec[1]}\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Keep only positions with variants in the given VCF file.")
	parser.add_argument("input_vcf", type = str, help = "Input VCF")
	parser.add_argument('--chr', type = str, required = False, default = None, help = "Chromosome ID")
	parser.add_argument('--bed', type = argparse.FileType("r"), required = False, default = None, help = "Region file")
	args = parser.parse_args()

	filter_positions(sys.stdin, args.input_vcf, args.bed, args.chr)
