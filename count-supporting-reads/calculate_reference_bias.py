#!/usr/bin/env python3

# Copyright (c) 2022-2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

import argparse
import sys


class ReadSupportParser(object):
	min_cov = 0
	is_first = True
	ref = ""
	alt = ""
	ref_count = 0
	alt_count = 0
	overall_ref_count = 0
	overall_alt_count = 0

	variants_counted = 0
	variants_skipped = 0


	def output_balance(self):
		if self.min_cov <= self.ref_count + self.alt_count:
			self.variants_counted += 1
			balance = float(self.ref_count) / float(self.ref_count + self.alt_count)
			ref_len = len(self.ref)
			alt_len = len(self.alt)
			sys.stdout.write(f"{balance}\t{ref_len}\t{alt_len}\n")
		else:
			self.variants_skipped += 1


	def parse(self, fp):
		sys.stdout.write(f"BALANCE\tREF_LENGTH\tALT_LENGTH\n")
		for lineno, line_ in enumerate(fp, start = 1):
			line = line_.rstrip("\n")
			fields = line.split("\t")

			rec_type = fields[0]
			if 'V' == rec_type:
				# Variant record.
				if not(self.is_first):
					self.output_balance()
				self.is_first = False

				alts = fields[5].split("\t")
				# FIXME: Handle this case properly? May not be needed b.c. we only have one diploid donor.
				if 1 != len(alts):
					sys.stderr.write(f"WARNING: ALT count is not equal to one on input line {lineno}. Considering only the first ALT.\n")
				self.alt = alts[0]
				if self.alt == "<DEL>":
					self.alt = ""
				self.ref = fields[4]

				is_reversed = int(fields[7])
				if is_reversed:
					self.ref, self.alt = self.alt, self.ref
				
				self.ref_count = 0
				self.alt_count = 0
			elif 'R' == rec_type:
				# Alignment (“read”) record.
				count = int(fields[1])
				text = fields[2]
				if text == self.ref:
					self.ref_count += count
					self.overall_ref_count += count
				elif text == self.alt:
					self.alt_count += count
					self.overall_alt_count += count
				else:
					#sys.stderr.write(f"WARNING: Allele “{text}” on line {lineno} matches neither ALT nor REF.\n")
					pass
			else:
				# Statistics, output as-is to stderr.
				sys.stderr.write(line_)

		if not(self.is_first):
			self.output_balance()

		sys.stderr.write(f"Variants counted: {self.variants_counted}\n")
		sys.stderr.write(f"Variants skipped: {self.variants_skipped}\n")
		if 0 < self.overall_alt_count:
			ratio = float(self.overall_ref_count) / float(self.overall_alt_count)
			sys.stdout.write(f"# Overall ref-to-alt ratio: {ratio}\n")
		else:
			sys.stderr.write("Found zero ALT alleles.\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'Calculate reference bias per site using #REF / (#REF + #ALT) from count_supporting_reads’s output.')
	parser.add_argument('--min-coverage', metavar = 'N', type = int, default = 1, help = 'minimum read coverage considering reads that support either REF or ALT')
	args = parser.parse_args()

	parser = ReadSupportParser()
	parser.min_cov = args.min_coverage
	parser.parse(sys.stdin)
