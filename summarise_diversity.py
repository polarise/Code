#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys

def main( fn ):
	with open( fn ) as f:
		genes = dict()
		for row in f:
#			print row.strip( "\n" )
			l = row.strip( "\n" ).split( "\t" )
#			print >> sys.stderr, l
			start = int( l[2].split( "(" )[0].split( "[" )[1] )
			if l[0] not in genes:
				genes[l[0]] = set([ start ])
			else:
				genes[l[0]].add( start )

	for gene,starts in genes.iteritems():
		print "\t".join([ gene, ",".join( map( str, list( starts ))), str( len( list( starts ))) ])

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <diversity_file>"
		sys.exit( 1 )
	
	
	main( fn )
