#!/usr/bin/env python
from __future__ import division
import sys
from BioClasses import *


def main( fn ):
	count_start_codons = dict()
	
	with open( fn ) as f:
		c = 0
		for row in f:
			if c > 100:
				break
			l = row.strip( "\n" ).split( "\t" )
			record = GTFRecord( *l )
			if record.feature == "start_codon":
				if record.group_dict['gene_id'] not in count_start_codons:
					count_start_codons[record.group_dict['gene_id']] = 1
				else:
					count_start_codons[record.group_dict['gene_id']] += 1
			c += 0
		
	for gene_id,count in count_start_codons.iteritems():
		if count > 1:
			print gene_id, count


if __name__ == "__main__":
	try:
		fn = sys.argv[1] # refseq gtf
	except IndexError:
		print >> sys.stderr, "usage:./script.py <gtf>"
		sys.exit( 0 )
	
	main( fn )
