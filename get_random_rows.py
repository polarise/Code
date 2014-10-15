#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import random

def main( fn, n ):
	with open( fn ) as f:
		all_lines = f.readlines()
		lines = [ random.choice( all_lines ) for i in xrange( n )]
	
	
	for l in lines:
		print l.strip( "\n" )

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <text-file> [<no-rows>]"
		sys.exit( 1 )
	
	try:
		n = int( sys.argv[2] )
	except IndexError:
		print >> sys.stderr, "Using default n=10..."
		n = 10
	
	main( fn, n )


