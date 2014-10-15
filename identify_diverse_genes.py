#!/usr/bin/env python
from __future__ import division
import sys
import pysam
from BioClasses import *
import cPickle
import argparse
import os

def main( fn, output_directory, genes_to_ignore, output_prefix, window_size, pvalue_thresh, excess, verbose, triplet_fix ):
	with open( fn ) as f:
		genes = cPickle.load( f )
	
	tabixfile = pysam.Tabixfile( "/home/paul/bioinf/Data/global_ribosome_profiles.bed.gz", "r" )
	
	peaks_file = open( "%s/%s.peaks.txt" % ( output_directory, output_prefix ), "w" )
	nopeaks_file = open( "%s/%s.nopeaks.txt" % ( output_directory, output_prefix ), "w" )
	
	c = 0
	d = 0
	for gene_id,G in genes.iteritems():
		if gene_id in genes_to_ignore:
			if verbose:
				print >> sys.stderr, "Ignoring gene: %s" % gene_id
			continue
		if c > 99:
			break
		equal_cds = G.get_equal_cds()
		if equal_cds is not None:
			for T in equal_cds:
				# get the regions for the
				# overal, UTR, CDS
				tx_region = T.region_str( True ) # True - zero-based
				utr_region = T.utr_region_str( True )
				cds_region = T.cds_region_str( True )
			
				tx_tabix_result = TabixResult( tx_region, tabixfile, verbose=verbose )
				UTR_tabix_result = TabixResult( utr_region, tabixfile, verbose=verbose )
				CDS_tabix_result = TabixResult( cds_region, tabixfile, verbose=verbose )
				
				# normalised to CDS including introns
#				UTR_tabix_result.compute_peaks( CDS_tabix_result.density, n=n, pvalue_thresh=pvalue_thresh, excess=excess )
				
				# normalised to CDS excluding introns
				UTR_tabix_result.compute_peaks( CDS_tabix_result.count/T.get_cds_length(), n=window_size, excess=excess, triplet_fix=triplet_fix )
				if len( UTR_tabix_result.found_peaks ) > 0:
					print >> peaks_file, "\t".join( map( str, [ G.gene_id, T, "%s|%s" % ( UTR_tabix_result.count, UTR_tabix_result.density ), "%s|%s" % ( CDS_tabix_result.count, round( CDS_tabix_result.count/T.get_cds_length(), 4 ) ), "%s|%s" % ( tx_tabix_result.count, tx_tabix_result.density ), ",".join([ "%s|%s" % ( pos, width ) for pos,width in UTR_tabix_result.found_peaks.iteritems()]) ]))
				else:
					print >> nopeaks_file, "\t".join( map( str, [ G.gene_id, T, "%s|%s" % ( UTR_tabix_result.count, UTR_tabix_result.density ), "%s|%s" % ( CDS_tabix_result.count, round( CDS_tabix_result.count/T.get_cds_length(), 4 ) ), "%s|%s" % ( tx_tabix_result.count, tx_tabix_result.density ), ",".join([ "%s|%s" % ( pos, width ) for pos,width in UTR_tabix_result.found_peaks.iteritems()]) ]))
		
				c += 0
		d += 1	

	tabixfile.close()
	peaks_file.close()
	nopeaks_file.close()
	
	print >> sys.stderr, "Results for %s of %s genes." % ( c, d )
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser( description="Identify diverse genes." )
	
	parser.add_argument( "-w", "--window-size", default=21, type=int, help="window size used to find peaks" )
	parser.add_argument( "-p", "--pvalue-thresh", default=0.01, type=float, help="p-value threshold to identify peaks" )
	parser.add_argument( "-e", "--excess-fold", default=5.0, type=float, help="amount by which mean of windowed region should exceed the CDS density" )
	parser.add_argument( "-v", "--verbose", action="store_true", default=False, help="print warnings to stderr" )
	parser.add_argument( "-i", "--ignore-genes", help="file containing genes to ignore" )
	parser.add_argument( "-o", "--output-prefix", default="results", help="prefix for output files" )
	parser.add_argument( "-d", "--output-directory", default=".", help="the directory into which to write the files; will be created if doesn't exist" )
	parser.add_argument( "-3", "--triplet-fix", default=False, action='store_true', help="this attempts to resolve peaks by accounting for triplet periodicity. It does this by considering the maximal count over a triplet window centered at the current nucleotide." )
	parser.add_argument( "file" )
	
	args = parser.parse_args()
	
	fn = args.file
	window_size = args.window_size
	pvalue_thresh = args.pvalue_thresh
	excess_fold = args.excess_fold
	verbose = args.verbose
	ignore_genes = args.ignore_genes
	output_prefix = args.output_prefix
	output_directory = args.output_directory
	triplet_fix = args.triplet_fix
	
	if window_size % 2 == 0:
		raise ValueError( "w (%s) must be odd!" % window_size )
	
	if ignore_genes is not None:
		with open( ignore_genes ) as f:
			if verbose:
				print >> sys.stderr, "Getting the list of files to ignore...",
			genes_to_ignore = { row.strip( "\n" ).split( "\t" )[0]:row.strip( "\n" ).split( "\t" )[1] for row in f }
			if verbose:
				print >> sys.stderr, "Done."
	else:
		genes_to_ignore = dict()
	
	if output_directory != ".":
		if verbose:
			print >> sys.stderr, "Attempting to create directory '%s%..." % output_directory
		try:
			os.mkdir( output_directory )
			if verbose:
				print >> sys.stderr, "Done."
		except OSError:
			print >> sys.stderr, "Directory exists... will overwrite existing files using same prefix (%s)." % output_prefix
	else:
		pass
	
	main( fn, output_directory, genes_to_ignore, output_prefix, window_size, pvalue_thresh, excess_fold, verbose, triplet_fix )
