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
			# tx_name	length	FS_bool
			# comp5317_c0_seq1	510	No
			l = row.strip( "\n" ).split( "\t" )
			stephen_txs[l[0]] = BioClasses.FrameshiftTranscript( name=l[0], \
					length=int( l[1] ))
	
	with open( path + f_all_FS_fn ) as f_all_FS:
		for row in f_all_FS:
			if row[0] == "T":
				continue
			# tx_name	length	#FS	FS_pos	FS_sig
			# comp7605_c0_seq1	5063	1	306	AAATAA
			l = row.strip( "\n" ).split( "\t" )
			if l[0] in stephen_txs:
				stephen_txs[l[0]].add_frameshift_site( int( l[3] ), l[4] )
	
	#for tx_id,T in stephen_txs.iteritems():
		#print T
	
	heuristic_txs = dict()
	f_heuristic_fn = "BioClasses/found_frameshifts.txt"
	with open( f_heuristic_fn ) as f_heuristic:
		for row in f_heuristic:
			l = row.strip( "\n" ).split( "\t" )
			# tx_name	length	loglik	FS_sig	f0	f1	FS_desig	FS_pos	pos_score	t0	t1	t2
			# comp1_c0_seq1	640	-816.906150722	AAATAG	0	1	+1	285	0.8906	0.4751	1.4394	2.1649
			if l[0] not in heuristic_txs:
				heuristic_txs[l[0]] = BioClasses.FrameshiftTranscript( name=l[0],\
					length=int( l[1] ))
				heuristic_txs[l[0]].add_frameshift_site( int( l[7] ), l[3] )
			else:
				heuristic_txs[l[0]].add_frameshift_site( int( l[7] ), l[3] )
	
	#for tx_id,T in heuristic_txs.iteritems():
		#print T
	
	print len( stephen_txs )
	print len( heuristic_txs )
	
	common = 0
	different = 0
	missing = 0
	for tx_id,T in stephen_txs.iteritems():
		try:
			if stephen_txs[tx_id].rough_equality( heuristic_txs[tx_id] ):
				common += 1
			else:
				different += 1
		except KeyError:
			missing += 1
	print common, different, missing
	print common + different, common + different + missing

if __name__ == "__main__":
	main()