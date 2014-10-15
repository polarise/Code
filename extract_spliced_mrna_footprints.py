#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
import pysam
import BioClasses

def main( bed_fn, gtf_fn ):
	# iterate through GTF and build gene models
	genes = dict()
	with open( gtf_fn ) as gtffile:
		c = 0
		for row in gtffile:
			if c > 1000:
				break
			l = row.strip( "\n" ).split( "\t" )
			record = BioClasses.GTFRecord( *l )
			
			if record.feature == "gene":
				genes[record.group_dict['gene_id']] = BioClasses.Gene( record )
			else:
				genes[record.group_dict['gene_id']].process_record( record )
			
			c += 0
	
	
	# tabix of BED
	tabixfile = pysam.Tabixfile( bed_fn, "r" )
	
	# for each gene
	transcripts = dict()
	for gene_id,G in genes.iteritems():
		# for each transcript
		for transcript_id,T in G.transcripts.iteritems():
			# get the data from the BED
			exon_numbers = T.exons.keys()
			exon_numbers.sort()
			for exon_number in exon_numbers:
				exon_region = T.exons[exon_number].region_str( zero_based=True)
				exon_tabix_result = BioClasses.TabixResult( exon_region, tabixfile )
				if transcript_id not in transcripts:
					if T.strand == "+":
						transcripts[transcript_id] = exon_tabix_result.serial_data
					elif T.strand == "-":
						transcripts[transcript_id] = exon_tabix_result.serial_data[::-1]
				else:
					if T.strand == "+":
						transcripts[transcript_id] += exon_tabix_result.serial_data
					elif T.strand == "-":
						transcripts[transcript_id] += exon_tabix_result.serial_data[::-1]
		
	for transcript_id,T in transcripts.iteritems():
		print "\t".join([transcript_id, ",".join( map( str, T ))])				

if __name__ == "__main__":
	gtf_fn = sys.argv[1]
	bed_fn = sys.argv[2]
	main( bed_fn, gtf_fn )
