# Copyright (c) 2023 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).

from Bio import SeqIO, SeqRecord
import argparse
import sys


def process(contig_ids, contig_id_output_fh, description):

	def helper():
		for rec in SeqIO.parse(sys.stdin, "fasta"):
			if contig_id_output_fh:
				contig_id_output_fh.write(rec.id)
				contig_id_output_fh.write('\n')
			
			if rec.id not in contig_ids:
				yield SeqRecord.SeqRecord(rec.seq, id = rec.id, description = description)
	
	def make_rec2title():
		if description is None:
			return lambda rec: rec.id
		else:
			return lambda rec: f"{rec.id}\t{rec.description}"
	
	writer = SeqIO.FastaIO.FastaTwoLineWriter(sys.stdout, make_rec2title())
	writer.write_file(helper())


if __name__ == "__main__":
	parser = argparse.ArgumentParser("Output reference sequences as read from stdin with the given entries removed.")
	parser.add_argument('-c', '--contig', type = str, nargs = '+', help = "Contig identifier to be removed")
	parser.add_argument('-d', '--description', type = str, help = "Record description")
	parser.add_argument('-o', '--output-contig-ids', type = argparse.FileType('w'), help = "Output contig identifiers from the reference")
	args = parser.parse_args()
	
	process(frozenset(args.contig), args.output_contig_ids, args.description)
