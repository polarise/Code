#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import pysam
import BioClasses

def PrintStatic( line ):
	sys.stderr.write( "\r%s".ljust( 20 ) % line )
	sys.stderr.flush()

def main( bed_fn, gtf_fn ):
	# iterate through GTF and build gene models
	genes = dict()
	print >> sys.stderr, "Reading GTF file...",
	with open( gtf_fn ) as gtffile:
		c = 0
		for row in gtffile:
			if c > 10000:
				break
			l = row.strip( "\n" ).split( "\t" )
			record = BioClasses.GTFRecord( *l )
			
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = BioClasses.Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )
			
			c += 0
	print >> sys.stderr, "Done."
	
	# tabix of BED
	print >> sys.stderr, "Opening BED file...",
	tabixfile = pysam.Tabixfile( bed_fn, "r" )
	print >> sys.stderr, "Done."
	
	no_of_genes = len( genes.keys() )
	
	# for each gene
	print >> sys.stderr, "Processing transcripts from %s genes..." % no_of_genes
	transcripts = dict()
	c = 0
	for gene_id,G in genes.iteritems():
		PrintStatic( "Progress: %.2f%%" % ( c/no_of_genes*100 ))
		#print >> sys.stderr, c, no_of_genes
		
		# for each transcript
		for transcript_id,T in G.transcripts.iteritems():
			# get the data from the BED
			#print >> sys.stderr, transcript_id+": ",
			exon_numbers = T.exons.keys()
			exon_numbers.sort()
			length = T.get_full_length()
			for exon_number in exon_numbers:
				#print >> sys.stderr, exon_number,
				exon_region = T.exons[exon_number].region_str( zero_based=True )
				exon_tabix_result = BioClasses.TabixResult( exon_region, tabixfile, verbose=False )
				if (transcript_id,T.strand,length) not in transcripts:
					transcripts[transcript_id,T.strand,length] = exon_tabix_result.serial_data
				else:
					if T.strand == "+":
						transcripts[transcript_id,T.strand,length] += exon_tabix_result.serial_data
					elif T.strand == "-":
						transcripts[transcript_id,T.strand,length] = exon_tabix_result.serial_data + transcripts[transcript_id,T.strand,length]
			#print >> sys.stderr
		c += 1
	PrintStatic( "Progress: %.2f%%" % ( c/no_of_genes*100 ))
	print >> sys.stderr
	
	for transcript_id,T in transcripts.iteritems():
		if transcript_id[1] == "+":
			print "\t".join([ transcript_id[0], transcript_id[1], str( transcript_id[2] ), str( len( T )), ",".join( map( str, T ))])
		elif transcript_id[1] == "-":
			print "\t".join([ transcript_id[0], transcript_id[1], str( transcript_id[2] ), str( len( T )), ",".join( map( str, T[::-1] ))])

if __name__ == "__main__":
	gtf_fn = sys.argv[1]
	bed_fn = sys.argv[2]
	main( bed_fn, gtf_fn )
