#! /usr/bin/env python


## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details

## Author(s): Kenzo-Hugo Hillion
## Contact(s): kenzo.hillion@curie.fr
## Python version: 2.7
## Script description: Script to select a read based on alignment score ("AS:i") mapped on both parental genomes

scriptVersion = '0.1 - 05-24-2016'

"""
Script to select a read based on alignment score ("AS:i") mapped on both parental genomes
	
	INPUT(s) :	
		BAM file of mapped reads sorted by name
	OUTPUT(s) :
		BAM file of selected reads with a new flag added :
		-XA:i:0 -> read is not assigned to any parental genome
		-XA:i:1 -> read is assigned to the first parental genome
		-XA:i:2 -> read is assigned to the second parental genome
		-XA:i:3 -> read is assigned to any parental genome because of mapping to different region with the same alignment score
	
	USE :
		./select_best.py -1 parent1.bam -2 parent2.bam -o output_dir/

"""

###########  Import  ###########

import sys
import getopt
import os
import re
import pysam
import itertools

########## Function(s) ##########

def usage():
	print "==== select_best.py. version " + scriptVersion + " ===="
 	print "Select best read between reads mapped on two parental genomes from the two BAM files sorted by names."

 	"""Usage function"""
 	print "Usage : compGen.py"
	print "-1/--par1 <names sorted BAM file of parent 1>"
	print "-2/--par2 <names sorted BAM file of parent 2>"
	print "[-f/--fused] <prefix name for output bam>"
	print "[-o/--output_dir] <Output directory>"
 	print "[-h/--help] <Help>"
	return

def get_args():
	"""Get argument"""
	try:
		opts, args = getopt.getopt(
 		sys.argv[1:],
 		"h:o:1:2:f:",
 		["help","output_dir","par1","par2","fused"])
	except getopt.GetoptError:
		usage()
		sys.exit(-1)
	return opts



########## Main ##########

if __name__ == "__main__":

	# Script arguments
	opts=get_args()
	if len(opts) == 0:
		usage()
		sys.exit()

	# Default values for arguments
	outname="selected"
	output_dir="./"
	par1=None
	par2=None

	# Store arguments
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-o", "--output_dir"):
			output_dir = arg
		elif opt in ("-1", "--par1"):
			par1 = arg
		elif opt in ("-2", "--par2"):
			par2 = arg
		elif opt in ("-f", "--fused"):
			outname = arg
		else:
			assert False, "unhandled option"

	# Test arguments
	if ((par1 is None) | (par2 is None)):
		sys.stderr.write("INPUT ERROR : BAM file(s) not specified in arguments.")
		usage()
		sys.exit(1)

	##################################################################

	# Read SAM/BAM files
	bampar1 = pysam.Samfile(par1, "rb")
	bampar2 = pysam.Samfile(par2, "rb")

	# Open a file for writing
	bamout = pysam.Samfile(output_dir + outname + ".bam", "wb", template=bampar1)

	# Set up counters and name of the tag
	counter_read = 0
	tag = "XA"

	# Comparing read one by one
	for readp1,readp2 in itertools.izip(bampar1.fetch(until_eof=True),bampar2.fetch(until_eof=True)):
		counter_read += 1					# Increments read counter
		
		if (readp1.flag == 4):				# If the read is unmapped, set the score to -1000 to make sure it looses the comparison to the score of a mapped read
			as_p1 = -1000
		else:
			as_p1 = readp1.get_tag("AS")
		if (readp2.flag == 4):				# If the read is unmapped, set the score to -1000 to make sure it looses the comparison to the score of a mapped read
			as_p2 =- 1000
		else:
			as_p2=readp2.get_tag("AS")
	
		if ((as_p1 == -1000) & (as_p2 == -1000)):	# If read is unmapped in both parent, set the flag to ambigous and write the first read
			readp1.set_tag(tag,0)
			bamout.write(readp1)
			continue

		if (readp1.pos == readp2.pos):		# If the reads mapped at the same position on the genome (Unlikely that the same read mapped at the same position in two different chromosome ?)
			if (as_p1 > as_p2):				# Best alignment score of parent 1: keep read mapped on p1
				#print ("1")
				readp1.set_tag(tag,1)
				bamout.write(readp1)
			elif (as_p1 < as_p2):			# Best alignment score of parent 2: keep read mapped on p2
				#print ("2")
				readp2.set_tag(tag,2)
				bamout.write(readp2)
			else: 							# Same alignment score: keep read mapped on p1 with flag "XA:i:0"
				#print ("3")
				readp1.set_tag(tag, 0)
				bamout.write(readp1)
		else: 								# Reads mapped at two different positions in each genome
			if (as_p1 > as_p2): 			# Best alignment score of parent 1: Keep the read mapped on p1
				#print ("4")
				readp1.set_tag(tag,1)
				bamout.write(readp1)
			elif (as_p1 < as_p2):			# Best alignment score of parent 2: Keep the read mapped on p2
				#print ("5")
				readp2.set_tag(tag,2)
				bamout.write(readp2)
			else: 							# Ambigous case, we set the flag to 3
				#print ("6")
				readp1.set_tag(tag,3)
				bamout.write(readp1)
		
		if (counter_read % 1000000 == 0):
			print ("Reads treated : " + str(counter_read))

	print ("Total number of reads : " + str(counter_read))

	# Closing SAM/BAM files
	bamout.close()
	bampar1.close()
	bampar2.close()
