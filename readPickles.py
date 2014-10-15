#!/usr/bin/env python
from __future__ import division
import cPickle
import sys,os,time
from re import search
import argparse

parser = argparse.ArgumentParser(description="Script to read PIC files")
parser.add_argument('fn',help="input PIC file")
parser.add_argument('-n','--no-items',type=int,help="number of lines")
parser.add_argument('-s','--sub-pic',help="make a sub-PIC file of the name given with the requested number of lines")
parser.add_argument('-c','--count',action='store_true',help="count the number of features in the PIC file")
parser.add_argument('-t','--type',action='store_true',default=False,help="show the types of key-values")

args = parser.parse_args()

fn = args.fn
no = args.no_items
sfn = args.sub_pic
count = args.count
the_type = args.type

if search(r".gz$",fn):
	import gzip
	f = gzip.open(fn)
else:
	f = open(fn)

data = cPickle.load(f)

if the_type:
	one_key = data.keys()[0]
	one_value = data.values()[0]
	print >> sys.stderr,"The types of the PIC are: Key is %s and Value is %s" % (type(one_key),type(one_value))
	if isinstance(one_value,list) or isinstance(one_value,tuple) or isinstance(one_value,dict):
		one_sub_value = one_value[0]
		print >> sys.stderr,"The Values have Subvalues of type %s" % type(one_sub_value)
	exit(0)

if count:
	print >> sys.stderr,"File contains %s items."%len(data.keys())
	f.close()
	exit(0)

if no:
	no_items = no
else:
	no = len(data.keys())+1
	no_items = "all"

print >> sys.stderr,"Printing %s items..."%no_items
print >> sys.stderr

c = 0
if sfn: sub_data = dict()
for d in data:
	if c > no: break
	print str(d)+"\t"+str(data[d])
	if sfn: sub_data[d] = data[d]
	c += 1

f.close()
if sfn: 
	f = open(sfn,'w')
	cPickle.dump(sub_data,f,cPickle.HIGHEST_PROTOCOL)
	f.close()
