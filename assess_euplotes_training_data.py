#!/usr/bin/env python
from __future__ import division
import sys
import cPickle
import BioClasses

def main():
	path = "/home/paul/Dropbox/Euplotes/FrameshiftPredictionData/"
	f_all_txs_fn = "All.Transcripts.csv"
	f_all_FS_fn = "Frameshift.Transcripts.csv"

	stephen_txs = dict()
	with open( path + f_all_txs_fn ) as f_all_txs:
		for row in f_all_txs:
			if row[0] == "T":
				continue
			# comp5317_c0_seq1	510	No
			l = row.strip( "\n" ).split( "\t" )
			stephen_txs[l[0]] = BioClasses.FrameshiftTranscript( name=l[0], \
					length=int( l[1] ))
	
	with open( path + f_all_FS_fn ) as f_all_FS:
		for row in f_all_FS:
			if row[0] == "T":
				continue
			# comp7605_c0_seq1	5063	1	306	AAATAA
			l = row.strip( "\n" ).split( "\t" )
			if l[0] in stephen_txs:
				stephen_txs[l[0]].add_frameshift_site( int( l[3] ), l[4] )
	
	for tx_id,T in stephen_txs.iteritems():
		print T

if __name__ == "__main__":
	main()