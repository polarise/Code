#!/usr/bin/env python
from __future__ import division
import sys
import BioClasses
import time

def myfunction( val ):
	return val+val

def main( fn ):
	with open( fn ) as f:
		lines = f.readlines()
	
	serial_start = time.time()
	result = map( myfunction, lines )
	serial_stop = time.time()
	
	parallel_time = dict()
	for p in xrange( 2, 25 ):
		parallel_start = time.time()
		PT = BioClasses.ParallelTask( lines, myfunction, num_of_processors=p )
		result = PT.run()
		parallel_stop = time.time()
		parallel_time[p] = parallel_stop - parallel_start
	
	print "Serial time:\t1-cpu\t%s" % ( serial_stop - serial_start )
	
	for p,t in parallel_time.iteritems():
		print "Parallel time:\t%s-cpus\t%s" % ( p, t )
	
	#for r in result:
		#print r

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <filename>"
		sys.exit( 1 )
		
	main( fn )
	