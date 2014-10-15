#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
import BioClasses
import cPickle

def main( fn ):
	with open( fn ) as f:
		genes = cPickle.load( f )
	
	c = 0
	for gene_id,G in genes.iteritems():
		if c > 100:
			break
		if G.is_protein_coding():
			for transcript_id,T in G.transcripts.iteritems():
				if T.is_complete():
					print gene_id,T,T.get_cds_length()
		
			c += 1

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <gtf-pic>"
		sys.exit( 1 )
	
	main( fn )
	
