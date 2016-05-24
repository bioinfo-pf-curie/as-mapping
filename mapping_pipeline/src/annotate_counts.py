#! /usr/bin/env python


## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details

## Author(s): Kenzo-Hugo Hillion, Eric Viara, Nicolas Servant
## Contact(s): kenzo.hillion@curie.fr
## Python version: 2.7
## Script description: Script to annotate the counts from check_variants (clinTools)

scriptVersion = '0.1 - 04-29-2016'

"""
This script is used to annotation the counts from check_variants.py (from clinTools)
It is based on a VCF file and associates the allele of the SNP to the correct strain
Some functions are from the markAllelicStatus.py script of HiC-Pro
    
    INPUT(s) : 
        Counts file from check_variants.py (clinTools) in the following format :
			barcode	chromosome	position	reference	depth	A	T	C	G	N	Ins	Del
		VCF file of differential SNPs between the two strains
    OUTPUT(s) :
        Table of counts per strain in the following format (ON STD_OUT) :
			BARCODE	CHR	POS	STRAIN2_BASE	STRAIN2_BASE	STRAIN1_COUNTS	STRAIN2_COUNTS	OTHER	TOTAL_COUNTS
    USE :
        ./annotate_counts.py -i counts.txt -s diff_snps.vcf > output.txt

"""
###########  Import  ###########

import getopt
import sys
import os
import re
import pysam
import time

###########  Function(s)  ###########

def usage():
	"""Usage function"""
	print "==== annotate_counts.py. version " + scriptVersion + " ===="
	print "Usage : annotate_counts.py"
	print "-i/--input <output file from check_variants.py (clinTools)>"
	print "-s/--snp <SNP file information - VCF format>"
	print "output table on stdout"
	print "[-h/--help] <Help>"
	return


def get_args():
	"""Get argument"""
	try:
		opts, args = getopt.getopt(
			sys.argv[1:],
			"i:s:h",
			["input=",
			 "snp=",
			 "help"])
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(-1)
	return opts


def get_snp_gt(gt, ref, alt):		  
	gtsnp = []
	
	## gtsnp.append(ref)
	snp_geno = re.split('/|\|', gt)
	## '.' are not considered
	if len(snp_geno) != 2:
		return [None,None]
	
	## First Allele
	if int(snp_geno[0]) == 0:
		gtsnp.append(ref)
	elif int(snp_geno[0]) == 1:
		gtsnp.append(alt)
	else:
		gtsnp.append(None)
	
	## Second Allele
	if int(snp_geno[1]) == 0:
		gtsnp.append(ref)
	elif int(snp_geno[1]) == 1:
		gtsnp.append(alt)
	else:
		gtsnp.append(None)

	return gtsnp
	


def load_vcf( in_file, filter_qual=False):
	"""
	Load a VCF file in a dict object
	
	in_file = path to VCF file [character]
	filter_qual = if True, only good quality SNV are selected [boolean]
	"""

	vcf_handle = open(in_file)	
	header = []
	samples = []
	snps = {}
	var_counter = 0
	snp_counter = 0
	for line in vcf_handle:
		line = line.rstrip()
		## for now we don't care about the header
		if line.startswith('##'):
			continue
		elif line.startswith('#'):
			header  = header = line.split('\t')			# Split by '\t' before, but the format of the generated vcf are with spaces
			header[0]	= header[0][1:]
			samples = [ s.split('.')[0] for s in header[9:] ]
			if len(samples) > 1:
				print >> sys.stderr, "Warning : Multisamples VCF detected. Only the first genotype will be used !"
			continue
		else:
			fields = line.split('\t',9)
			var_counter+=1
			n = len(fields)
			chrom = "chr" + fields[0]
			start = int(fields[1])
			ref = fields[3]
			alt = fields[4]
			qfilter = fields[6]
			## Check format for first variant
			if var_counter == 1:
				format = fields[8] if n>8 else None
				if format.split(':')[0] != "GT":
					print >> sys.stderr,"Error : Invalid format - GT not detected at first position in ", format		 
					sys.exit(-1)

			genotypes  = fields[9].split('\t') if fields[9] else []
			geno = get_snp_gt(genotypes[0].split(':')[0], ref, alt)	  
			if filter_qual == False or (filter_qual == True and qfilter=="PASS"):
				## store only discriminant SNP
				if geno[0] != geno[1]:
					snp_counter+=1
					snps[(str(chrom), int(start), '1')] = geno[0]
					snps[(str(chrom), int(start), '2')] = geno[1]

	vcf_handle.close()
	return snps

def get_strains(in_file):
	"""
	Get the name of the two strains from the VCF file

	in_file = path to VCF file [character]
	"""
	strains=[]
	header=[]
	vcf_handle = open(in_file)
	for line in vcf_handle:
		line = line.rstrip()
		if line.startswith('##'):
			continue
		elif line.startswith('#'):
			header = line.split('\t')
			strains = header[9].split('-')
			break
		else:
			break
	vcf_handle.close()
	return strains
	
###########  Main  ###########

if __name__ == "__main__":

	# Script arguments
	opts = get_args()
	if len(opts) == 0:
		usage()
		sys.exit()
	
	# Default values for arguments
	inputFile = None
	snpFile = None
	
	# Store arguments
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-s", "--snp"):
			snpFile = arg
		elif opt in ("-i", "--ibam"):
			inputFile = arg
		else:
			assert False, "unhandled option"

	# Test arguments
	if (snpFile is None):
        print "INPUT ERROR : VCF file not specified in arguments."
        usage()
        sys.exit()
	if (inputFile is None):
        print "INPUT ERROR : Counts file not specified in arguments."
        usage()
        sys.exit()

	##################################################################
	
	# Read the SNP file
	snps = {}
	snps = load_vcf(snpFile, filter_qual=False)

	# Variables initialization
	strains = []
	barcode = ""
	chromosome = ""
	position = 0
	geno1 = ""	
	geno2 = ""
	count1 = 0
	count2 = 0
	other = 0
	total = 0
	headpos = {}			# The dictionnary contains the order of every element of the input file

	# Print header of the file on standard output
	strains=get_strains(snpFile) # strains[0] corresponds to the first genotype and strains[1] to the second
	print "BARCODE\tCHR\tPOS\t" + strains[0] + "_BASE\t" + strains[1] + "_BASE\t" + strains[0] \
 + "_COUNTS\t" + strains[1] + "_COUNTS\tOTHER\tTOTAL"
	
	
	file_handle = open(inputFile)
	for line in file_handle:
		line = line.rstrip()
		
		# if the line corresponds to the header of the count file
		if line.startswith('barcode'):
			header = line.split('\t')
			cpt=0
			for element in header:
				headpos[str(element)]=cpt
				cpt+=1
			continue	
		
		# We now process the other lines of the file (variables are used to make the code clear)
		fields = line.split('\t')
		barcode=str(fields[0])
		chromosome=str(fields[1])
		position=int(fields[2])
		# We get the genotypes of the strains from the dictionnary
		geno1=snps[(chromosome,position,'1')]
		geno2=snps[(chromosome,position,'2')]
		# We look at the counts for the corresponding allele
		count1=int(fields[headpos[geno1]])
		count2=int(fields[headpos[geno2]])
		total=int(fields[4])
		other=total-(count1+count2)
		# Write the line on standard output
		print barcode + "\t" + chromosome + "\t" + str(position) + "\t" + geno1 + "\t" + geno2 + "\t" + str(count1) + "\t" + str(count2) + "\t" + str(other) + "\t" + str(total)

	# Close handler
	file_handle.close()
