#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import pysam
import BioClasses
import cPickle

def length_from_region( region ):
	seqname, start_stop = region.split( ":" )
	start, stop = map( int, start_stop.split( "-" ))
	return stop - start + 1

def main( fn ):
	with open( fn ) as f:
		genes = cPickle.load( f )
	
	print "\t".join([ "gene_id", "UTR_len", "CDS_len", "CDS_len_intronless" ])
	for gene_id,G in genes.iteritems():
		equal_cds = G.get_equal_cds()
		if equal_cds is not None:
			print "\t".join( map( str, [ gene_id, length_from_region( G.get_overall_UTR_region( "equal_cds" )), length_from_region( equal_cds[0].cds_region_str()), equal_cds[0].get_cds_length() ]))

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <pic-file>"
		sys.exit( 1 )
	main( fn )
