#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import pysam
import BioClasses

def main():
	tabixfile = pysam.Tabixfile( "/home/paul/bioinf/Data/global_ribosome_profiles.1-based.bed.gz", "r" )
	
	region = "chr1:1187065-1187832"
	region = "chr1:1190849-1190868"
	region = "chr1:1190819-1190868"
	region = "chr1:11869-12227"
	region = "chr1:93299173-93299203"
	
	tabix_result = BioClasses.TabixResult( region, tabixfile )
	
	print region
	print tabix_result.serial_data
	print tabix_result.length

if __name__ == "__main__":
	main()