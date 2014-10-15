#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
import gzip
import pysam
import BioClasses
import cPickle

def main( gtf_fn, genes_fn ):
	# genes
	genes = dict()

	with open( gtf_fn ) as f:
		c = 0
		for row in f:
			if c > 10000:
				break
			
			l = row.strip( "\n" ).split( "\t" )
			record = BioClasses.GTFRecord( *l )
		
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = BioClasses.Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )
		
			c += 0
	
#	for transcript_id,T in genes['ENSG00000211769'].transcripts.iteritems():
#		print transcript_id,T.is_complete()
#	
#	sys.exit( 0 )	
	
	# get the tabix file handler
	tabixfile = pysam.Tabixfile( genes_fn, "r" )
	
	# iterate over the genes
	d = 0
	for gene_id,G in genes.iteritems():
		if d > 1000:
			break
		# get the overall UTR for equal_cds transcripts for this gene
		utr_region = G.get_overall_UTR_region( "equal_cds" )
		
		if utr_region is None:
			print >> sys.stderr, "No region for %s" % gene_id
			continue
		
		# get overlaps
		tabix_result = BioClasses.TabixResult( utr_region, tabixfile, filetype="gtf" )
		
		# if there are overlaps...
		if tabix_result.record_count > 1:
			for record in tabix_result.records:
				print "\t".join([ gene_id, utr_region, record.group_dict['gene_id'], record.region_str() ]) # may result in duplicates
		d += 0
					
	tabixfile.close()


if __name__ == "__main__":
	try:
		gtf_fn = sys.argv[1]
		genes_fn = sys.argv[2]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <full-gtf> <bgzipped and indexed gene.gtf>"
		sys.exit( 1 )
	
	main( gtf_fn, genes_fn )
