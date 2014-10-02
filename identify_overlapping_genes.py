#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
import gzip
import pysam
import BioClasses

def main( gene_fn ):
	# get the tabix file handler
	tabixfile = pysam.Tabixfile( gene_fn, "r" )
	
	# iterate over the bgzip GTF file
	with gzip.open( gene_fn ) as f:
		for row in f:
			l = row.strip( "\n" ).split( "\t" )
			record = BioClasses.GTFRecord( *l )
			region = "%s:%s-%s" % ( record.seqname, record.start, record.end )
			tabix_result = BioClasses.TabixResult( region, tabixfile, filetype="gtf" )
			if tabix_result.record_count == 1:
				for record in tabix_result.records:
					print "\t".join([ record.group_dict['gene_id'], record.region_str() ])
			
	tabixfile.close()


if __name__ == "__main__":
	try:
		genes_fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <bgzipped and indexed gene.gtf>"
		sys.exit( 1 )
	
	main( genes_fn )
