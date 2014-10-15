#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
from BioClasses import *

def main( fn ):
	genes = dict()
	
	with open( fn ) as f:
		c = 0
		for row in f:
			if c > 20000:
				break
				
			l = row.strip( "\n" ).split( "\t" )
			record = GTFRecord( *l )
			
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )
			
			c += 1
		
	for gene_id,G in genes.iteritems():
		equal_cds = G.get_equal_cds()
		if len( equal_cds ) > 1:
			print "\t".join([ gene_id, G.region_str(), G.strand, G.get_overall_UTR_region( "equal_cds" ), "|".join( map( lambda x: x.transcript_id, equal_cds ))])

if __name__ == "__main__":
	try:
		fn = sys.argv[1] # gtf-file
	except IndexError:
		print >> sys.stderr, "usage:./script.py <gtf-file>"
		sys.exit( 1 )
	
	main( fn )
