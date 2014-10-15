#!/usr/bin/env python
from __future__ import division
import sys
import BioClasses
import cPickle

def main( fn ):
	genes = dict()
	with open( fn ) as f:
		for row in f:
			l = row.strip( "\n" ).split( "\t" )
			record = BioClasses.GTFRecord( *l )
			
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = BioClasses.Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )

	with open( "/home/paul/bioinf/Resources/H_sapiens/hg19.chr.pic", "w" ) as f:	
		cPickle.dump( genes, f, cPickle.HIGHEST_PROTOCOL )

if __name__ == "__main__":
	try:
		fn = sys.argv[1] # refseq gtf
	except IndexError:
		print >> sys.stderr, "usage:./script.py <gtf>"
		sys.exit( 0 )
	
	main( fn )
