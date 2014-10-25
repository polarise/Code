#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import pysam
import BioClasses
from BioClasses.Utils import *
import cPickle
import argparse
import os
import time
import multiprocessing
		
def get_peaks( gene, tabixfile, lock ):
	out_str = str()
	equal_cds = gene.get_equal_cds()
	if equal_cds is not None:
		for T in equal_cds:
			# get the regions for the
			# overal, UTR, CDS
			tx_region = T.region_str( True ) # True - zero-based
			utr_region = T.utr_region_str( True )
			cds_region = T.cds_region_str( True )
			
			lock.acquire()
			tx_tabix_result = BioClasses.TabixResult( tx_region, tabixfile, verbose=False )
			UTR_tabix_result = BioClasses.TabixResult( utr_region, tabixfile, verbose=False )
			CDS_tabix_result = BioClasses.TabixResult( cds_region, tabixfile, verbose=False )
			lock.release()
			
			# normalised to CDS excluding introns
			UTR_tabix_result.compute_peaks( CDS_tabix_result.count/T.get_cds_length(), n=9, excess=2, triplet_fix=True )
			if len( UTR_tabix_result.found_peaks ) > 0:
				out_str = "\t".join( map( str, [ "peak", gene.gene_id, T, "%s|%s" % ( UTR_tabix_result.count, UTR_tabix_result.density ), "%s|%s" % ( CDS_tabix_result.count, round( CDS_tabix_result.count/T.get_cds_length(), 4 ) ), "%s|%s" % ( tx_tabix_result.count, tx_tabix_result.density ), ",".join([ "%s|%s" % ( pos, width ) for pos,width in UTR_tabix_result.found_peaks.iteritems()]) ]))
			else:
				out_str = "\t".join( map( str, [ "nopeak", gene.gene_id, T, "%s|%s" % ( UTR_tabix_result.count, UTR_tabix_result.density ), "%s|%s" % ( CDS_tabix_result.count, round( CDS_tabix_result.count/T.get_cds_length(), 4 ) ), "%s|%s" % ( tx_tabix_result.count, tx_tabix_result.density ), ",".join([ "%s|%s" % ( pos, width ) for pos,width in UTR_tabix_result.found_peaks.iteritems()]) ]))
	return out_str

def main( fn ):
	msg( "Loading data from %s..." % fn, False ) 
	with open( fn ) as f:
		genes = cPickle.load( f )
	msg( "Done." )
	msg( "Found gene objects for %s genes." % len( genes ))
	
	bed_fn = "/home/paul/bioinf/Data/global_ribosome_profiles.1-based.bed.gz"
	msg( "Opening tabix file %s..." % bed_fn, False )
	tabixfile = pysam.Tabixfile( bed_fn, "r" )
	msg( "Done." )
	
	msg( "Starting parallel task on 20 CPUS...", False )
	L = multiprocessing.Lock()
	PT = BioClasses.ParallelTask( genes.values(), callback=get_peaks, num_of_processors=20, var=tabixfile, lock=L )
	results = PT.run()
	msg( "Done." )
	
	msg( "Writing to files...", False )
	peaks_file = open( "%s/%s.peaks.txt" % ( output_directory, output_prefix ), "w" )
	nopeaks_file = open( "%s/%s.nopeaks.txt" % ( output_directory, output_prefix ), "w" )
	
	for r in results:
		l = r.split( "\t" )
		if l[0] == "peak":
			print >> peaks_file, "\t".join( l[1:] )
		elif l[0] == "nopeak":
			print >> nopeaks_file, "\t".join( l[1:] )
	
	peaks_file.close()
	nopeaks_file.close()
	
	msg( "Done." )
	msg( "Results for %s of %s genes." % ( c, d ))
	
if __name__ == "__main__":
	main( sys.argv[1] ) # sys.argv[1] is the GTF file