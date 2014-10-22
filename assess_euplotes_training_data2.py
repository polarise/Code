#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import scipy
import argparse
import BioClasses

class StephenFrameshifts( object ):
	def __init__( self, filename ):
		self.transcripts  = dict()
		with open( filename ) as f:
			for row in f:
				if row[0] == "T":
					continue
				l = row.strip( "\n" ).split( "\t" )
				self.transcripts[l[0]] = BioClasses.FrameshiftTranscript( name=l[0], \
					length=int( l[1] ))
	
	def add_true_frameshifts( self, filename ):
		with open( filename ) as f:
			for row in f:
				if row[0] == "T":
					continue
				# tx_name	length	#FS	FS_pos	FS_sig
				# comp7605_c0_seq1	5063	1	306	AAATAA
				l = row.strip( "\n" ).split( "\t" )
				if l[0] in self.transcripts:
					self.transcripts[l[0]].add_frameshift_site( int( l[3] ), l[4] )
	
	def __getitem__( self, item ):
		try:
			return self.transcripts[item]
		except KeyError:
			return None

class HeuristicFrameshifts( object ):
	def __init__( self, filename ):
		self.transcripts = dict()
		with open( filename ) as f:
			for row in f:
				l = row.strip( "\n" ).split( "\t" )
				# tx_name	length	loglik	FS_sig	f0	f1	FS_desig	FS_pos	pos_score	t0	t1	t2
				# comp1_c0_seq1	640	-816.906150722	AAATAG	0	1	+1	285	0.8906	0.4751	1.4394	2.1649
				if l[0] not in self.transcripts:
					self.transcripts[l[0]] = BioClasses.FrameshiftTranscript( name=l[0],\
						length=int( l[1] ))
					self.transcripts[l[0]].add_frameshift_site( int( l[7] ), l[3], ( float( l[9] ), float( l[10] ), float( l[11] )))
				else:
					self.transcripts[l[0]].add_frameshift_site( int( l[7] ), l[3], ( float( l[9] ), float( l[10] ), float( l[11] )))
		
	def add_single_frames( self, filename ):
		with open( filename ) as f:
			for row in f:
				l = row.strip( "\n" ).split( "\t" )
				self.transcripts[l[0]] = BioClasses.FrameshiftTranscript( name=l[0],\
					length=1 )
		
	def __getitem__( self, item ):
		try:
			return self.transcripts[item]
		except KeyError:
			return None
	
	def simple_count( self, p0, p1, theta0, refFST=None ):
		"""
		p0 - starting position in %
		p1 - final position in %
		"""
		TP = 0
		TN = 0
		FP = 0
		FN = 0
		if refFST != None: # we're checking against some reference
			for tx,FST in refFST.transcripts.iteritems():
				heuristicFST = self[tx]
				if heuristicFST != None: # it is present in heuristics
					if heuristicFST.has_frameshift( p0, p1, theta0 ) and FST.has_frameshift( p0, p1 ):
						TP += 1
					elif heuristicFST.has_frameshift( p0, p1, theta0 ) and not FST.has_frameshift( p0, p1 ):
						FP += 1
					elif not heuristicFST.has_frameshift( p0, p1, theta0 ) and FST.has_frameshift( p0, p1 ):
						FN += 1
					elif not heuristicFST.has_frameshift( p0, p1, theta0 ) and not FST.has_frameshift( p0, p1 ):
						TN += 1
			return TP,TN,FP,FN
		else:
			pass
	
	def complex_count( self, p0, p1, theta0, refFST=None, tol=3 ):
		TP = 0
		TN = 0
		FP = 0
		FN = 0
		if refFST != None:
			for tx,FST in refFST.transcripts.iteritems():
				heuristicFST = self[tx]
				if heuristicFST != None:
					if heuristicFST.has_frameshift( p0, p1, theta0 ) and FST.has_frameshift( p0, p1 ):
						if heuristicFST.has_exact_frameshift( FST, p0, p1, theta0, tol=tol ):
							#print "Heuristic:"
							#print heuristicFST.filtered_print( p0, p1, theta0 )
							#print "Stephen:"
							#print FST
							#print "----------------------------------------"
							TP += 1
						else:
							FP += 1
					elif heuristicFST.has_frameshift( p0, p1, theta0 ) and not FST.has_frameshift( p0, p1 ):
						FP += 1
					elif not heuristicFST.has_frameshift( p0, p1, theta0 ) and FST.has_frameshift( p0, p1 ):
						FN += 1
					elif not heuristicFST.has_frameshift( p0, p1, theta0 ) and not FST.has_frameshift( p0, p1 ):
						TN += 1
			return TP,TN,FP,FN
		else:
			pass					
	
	def filter_frameshifts( self, p0, p1, theta0 ):
		fss_list = list()
		for fst_i,fst in self.transcripts.iteritems():
			for fss in fst.frameshifts( p0, p1, theta0 ):
				fss_list.append( fss )
		
		return fss_list

#===============================================================================

def main( begin, end, simple_count, complex_count ):
	sfs = StephenFrameshifts( "/home/paul/Dropbox/Euplotes/FrameshiftPredictionData/All.Transcripts.csv" )
	sfs.add_true_frameshifts( "/home/paul/Dropbox/Euplotes/FrameshiftPredictionData/Frameshift.Transcripts.csv" )
	hfs = HeuristicFrameshifts( "BioClasses/found_frameshifts.txt" )
	hfs.add_single_frames( "BioClasses/single_frame.txt" )

	for t0 in scipy.linspace( 0, 1.9, 101 ): # main parameter that selects FS
		#t0 = 0.911061869541
		# we use D1 as a reference
		if simple_count:
			TP, TN, FP, FN = hfs.simple_count( refFST=sfs, p0=begin, p1=end, theta0=t0 )
		elif complex_count:
			TP, TN, FP, FN = hfs.complex_count( refFST=sfs, p0=begin, p1=end, theta0=t0, tol=2 )
		TPR = TP/( TP + FN ) # true positive rate
		FPR = FP/( TN + FP ) # false positive rate
		ACC = ( TP + TN )/( TP + TN + FP + FN )	# accuracy
		try:
			PPV = TP/( TP + FP ) # precision (positive predictive value)
		except ZeroDivisionError:
			PPV = None
		print "\t".join( map( str, [ begin, end, t0, TP, FP, TN, FN, TP+TN+FP+FN, TPR, FPR, ACC, PPV ]))
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser( description="Script to compute TPR, FPR and ACC." )
	
	parser.add_argument( "-0", "--begin", default=0, type=float, help="normalised beginning of transcript [default: 0]" )
	parser.add_argument( "-1", "--end", default=1, type=float, help="normalised end of transcript [default: 1]" )
	parser.add_argument( "-s", "--simple-count", default=False, action='store_true', help="simple count i.e. TP is when there is detection without detail [default:True]" )
	parser.add_argument( "-c", "--complex-count", default=False, action='store_true', help="simple count i.e. TP is when there is detection with detail [default:False]" )
	
	args = parser.parse_args()
	
	begin = args.begin
	end = args.end
	try:
		assert 0 <= begin < end <= 1
	except:
		raise ValueError( "Invalid values for 'begin' (%s) and/or 'end' (%s)." % ( begin, end ))
	
	simple_count = args.simple_count
	complex_count = args.complex_count
	
	if simple_count and complex_count:
		print >> sys.stderr, "You must specify only one of --simple-count or --complex-count."
		print >> sys.stderr, "Running with --simple-count by default."
		complex_count = False
	elif not simple_count and not complex_count:
		simple_count = True
	
	main( begin, end, simple_count, complex_count )
		