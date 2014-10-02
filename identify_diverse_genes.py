#!/usr/bin/env python
from __future__ import division
import sys
import pysam
from BioClasses import *


def main( fn, n ):
	genes = dict()
	
	with open( fn ) as f:
		c = 0
		for row in f:
			if c > 10000:
				break
				
			l = row.strip( "\n" ).split( "\t" )
			record = GTFRecord( *l )
			
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )
			
			c += 0
	
	tabixfile = pysam.Tabixfile( "/home/paul/bioinf/Data/global_ribosome_profiles.bed.gz", "r" )	
	
#	for gene_id,G in genes.iteritems():
#		protein_coding = G.get_protein_coding()
#		equal_cds = G.get_equal_cds()
#		if len( equal_cds ) > 1:
#			print "%s\t%s\t%s\t%s" % ( G, str( len( equal_cds )), str( len( protein_coding )), ",".join( map( str,  equal_cds )))
	
	for gene_id,G in genes.iteritems():
		equal_cds = G.get_equal_cds()
		if len( equal_cds ) > 1:
			for T in equal_cds:
				# get the regions for the
				# overal, UTR, CDS
				tx_region = T.region_str()
				utr_region = T.utr_region_str()
				cds_region = T.cds_region_str()
			
				tx_tabix_result = TabixResult( tx_region, tabixfile )
				UTR_tabix_result = TabixResult( utr_region, tabixfile )
				CDS_tabix_result = TabixResult( cds_region, tabixfile )
				
				UTR_tabix_result.compute_peaks( CDS_tabix_result.density, n=n )

				print "\t".join( map( str, [ G.gene_id, T, "%s|%s" % ( UTR_tabix_result.count, UTR_tabix_result.density ), "%s|%s" % ( CDS_tabix_result.count, CDS_tabix_result.density ), "%s|%s" % ( tx_tabix_result.count, tx_tabix_result.density ), ",".join([ "%s|%s" % ( pos, width ) for pos,width in UTR_tabix_result.found_peaks.iteritems()]) ]))				
				
#				print tx_tabix_result.serial_data, len( tx_tabix_result.serial_data )
#				print sum( UTR_tabix_result.serial_data ), len( UTR_tabix_result.serial_data )
#				print CDS_tabix_result.serial_data, len( CDS_tabix_result.serial_data )
#				print sum( UTR_tabix_result.peak_data ), len( UTR_tabix_result.peak_data )
#				print sum( UTR_tabix_result.peak_pvalues ), len( UTR_tabix_result.peak_pvalues )
#				print ",".join([ "%s|%s" for pos,width UTR_tabix_result.found_peaks.iteritems()])
#				print

	tabixfile.close()
	
if __name__ == "__main__":
	try:
		fn = sys.argv[1] # refseq gtf
	except IndexError:
		print >> sys.stderr, "usage:./script.py <gtf>"
		sys.exit( 0 )
	
	try:
		n = int( sys.argv[2] )
	except IndexError:
		print >> sys.stderr, "Warning: using default n of 21"
		n = 21
	
	main( fn, n )
