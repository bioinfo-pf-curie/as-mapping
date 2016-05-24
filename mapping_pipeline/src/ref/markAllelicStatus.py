#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie								
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## --- Modification by Kenzo-Hugo Hillion // 05-24-2016 ---


"""
Kenzo-Hugo Hillion 05-24-2016 :
The script has been adapted to mark the Allelic status
of reads mapped on the reference genome based on the presence
of SNPs in every reads.
The script goes through all the bases of every read and check if the
position is contained in the dictionnary.

Function added :
	-get_snps_positions(chrom, read, snps);

Function_removed :
	-get_mismatches_positions(read, base=None);
	-getGenomePos(read, pos);
	-get_read_tag(read, tag);
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
	print "Usage : python markAllelicStatus.py"
	print "-i/--ibam <BAM/SAM file of mapped reads>"
	print "-s/--snp <SNP file information - VCF format>"
	print "[-r/--rstat] <Generate a report with descriptive statistics>"
	print "[-t/--tag] <tag used to report allelic expression. Default XA>"
	print "[-o/--out] <output BAM file. Default is stdin>"
	print "[-v/--verbose] <verbose>"
	print "[-h/--help] <Help>"
	return


def get_args():
	"""Get argument"""
	try:
		opts, args = getopt.getopt(
			sys.argv[1:],
			"i:s:t:o:rvh",
			["ibam=",
			 "snp=",
			 "tag=",
			 "output=",
			 "rstat", "verbose", "help"])
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
	


def load_vcf( in_file, filter_qual=False, verbose=False, debug=False ):
	"""
	Load a VCF file in a dict object
	
	in_file = path to VCF file [character]
	filter_qual = if True, only good quality SNV are selected [boolean]
	verbose = if True, verbose mode[boolean]
	debug = if True, debug mode [boolean]
	"""
	if verbose:
		print "## Loading VCF file '", in_file, "'..."

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
			chrom = fields[0]
			start = int(fields[1])-1 ## 0-based
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
				if debug:
					print >> sys.stderr, str(chrom) + " - " + str(start) + " - "+ str(qfilter) +" -REF= " + str(ref) + " -ALT= " + str(alt) + " - G1=" + str(geno[0]) + " - G2=" + str(geno[1])
				## store only discriminant SNP
				if geno[0] != geno[1]:
					snp_counter+=1
					snps[(str(chrom), int(start), '1')] = geno[0]
					snps[(str(chrom), int(start), '2')] = geno[1]

		if (var_counter % 100000 == 0 and verbose):
				print "##", var_counter
 
	vcf_handle.close()
	if verbose:
		print "## Number of loaded SNPs =",len(snps)/2, "over",var_counter
	return snps


def getBaseAt(read, pos):
	"""
	Extract nucleotide within a read at a given set of positions

	read = read object [class pysam.AlignedSegment]
	pos = positions to convert [list]

	"""
	nuc = []
	for p in pos:
		#print (p)
		nuc.append(read.seq[(p-read.pos)])
	return nuc


def getAllelicStatus(chrom, gpos, genotype, snps, debug=False):
	"""
	For a given set of genomic position and assoctiated genotype, compare to a snp file and return a code status
	0 : unassigned - no snp information extracted from the read
	1 : genotype from REF genome is found
	2 : genotype from ALT genome is found
	3 : conflicting information extracted from the read

	gpos = genomic position to look at [list]
	genotype = read genotype [list]
	snps = snp information extracted from VCF file
	"""

	code = None
	g1_count = 0
	g2_count = 0
	l = len(genotype)
	chrn = re.sub("^[Cc]hr","",chrom)

	for i in range(len(genotype)):
		#print >> sys.stderr, chrn, gpos[i], genotype[i]
		if gpos[i] != None:
			if snps.has_key((str(chrn), int(gpos[i]), '1')) and snps.has_key((str(chrn), int(gpos[i]), '2')):
				if snps[(str(chrn), int(gpos[i]), '1')] == genotype[i]:
					g1_count+=1
				elif snps[(str(chrn), int(gpos[i]), '2')] == genotype[i]:
					g2_count+=1
				else:
					print >> sys.stderr, "Warning : no SNPs found at position " + chrom + ":" + str(gpos[i]+1) + ". SNP ignored"

	if g1_count > 0 and g2_count > 0:
		code = 3
	elif g1_count > 0 and g2_count == 0:
		code = 1
	elif g2_count > 0 and g1_count == 0:
		code = 2
	elif g1_count == 0 and g2_count == 0:
		code = 0

	return code


def get_snps_positions(chrom, read, snps):
	"""
	Extract the snp(s) positions within the read
	
	read = read object [class pysam.AlignedSegment]
	snps = snp information extracted from VCF file

	"""
	chrn = re.sub("^[Cc]hr","",chrom)
	genomePos = read.pos
	npos = [] ## 0-based
	# For every position of the read, keep the position if the key exists in our SNP dictionnary
	for pos in range(genomePos,genomePos+100):
		if snps.has_key((str(chrn), pos, '1')) and snps.has_key((str(chrn), pos, '2')):
			npos.append(pos)
	return npos


###########  Main  ###########

if __name__ == "__main__":


	# Read command line arguments
	opts = get_args()
	inputFile = None
	verbose = False
	report = False
	debug = False
	output = "-"
	tag = "XA"
	snps={}
	
	if len(opts) == 0:
		usage()
		sys.exit()

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-s", "--snp"):
			snpFile = arg
		elif opt in ("-i", "--ibam"):
			mappedReadsFile = arg
		elif opt in ("-t", "--tag"):
			tag = arg
		elif opt in ("-o", "--out"):
			output = arg
		elif opt in ("-r", "--rstat"):
			report = True
		elif opt in ("-v", "--verbose"):
			verbose = True
		else:
			assert False, "unhandled option"

	# Initialize variables
	reads_counter = 0
	ua_counter = 0
	g1_counter = 0
	g2_counter = 0
	cf_counter = 0
	SNP_counter = 0

	# Read the SNP file
	snps = load_vcf(snpFile, filter_qual=False, verbose=verbose, debug=False)
	
	# Read the SAM/BAM file
	if verbose:
		print "## Opening SAM/BAM file '", mappedReadsFile, "'..."
	infile = pysam.Samfile(mappedReadsFile, "rb")

	#samOut:
	if output == "-":
		outfile = pysam.AlignmentFile(output, "w", template=infile)
	else:
		outfile = pysam.AlignmentFile(output, "wb", template=infile)
	 
	# Verbose mode																																						
	if verbose:
		print "## " + __file__
		print "## ibam=", mappedReadsFile
		print "## snpFile=", snpFile
		print "## tag=", tag
		print "## output=" + output 
		print "## verbose=", verbose, "\n"

	# Reads are 0-based too (for both SAM and BAM format)
	# Loop on all reads
	if verbose:
		print "## Assigning allele specific information ..."

	for read in infile.fetch(until_eof=True):
		reads_counter += 1
		if not read.is_unmapped:## and read.cigarstring.find("D") != -1:
			read_chrom = infile.getrname(read.tid)
			SNPreadpos = get_snps_positions(read_chrom, read, snps) 	# Get the positions of the SNPs where the read mapped
			if (len(SNPreadpos)>0):
				#SNP_counter += len(Nreadpos)
				SNP_counter += 1 		# Modified instead of previous line
				SNPbase = getBaseAt(read, SNPreadpos)
				tagval = getAllelicStatus(read_chrom, SNPreadpos, SNPbase, snps, debug=debug)
				read.set_tag(tag, tagval)
			
				if debug:
					for i in range(len(SNPreadpos)):
						if Ngenomepos[i] != None:
							print >> sys.stderr, str(read_chrom) +"\t"+ str(Ngenomepos[i]) + "\t" + str(Ngenomepos[i]+1) + "\t" + str(read.qname) + "/N/" + str(SNPbase[i]) + "\t" + str(tagval)
				if tagval == 0:
					ua_counter += 1
				elif tagval == 1:
					g1_counter += 1
				elif tagval == 2:
					g2_counter += 1
				elif tagval == 3:
					cf_counter += 1										
			else:
				read.set_tag(tag, 0)
				ua_counter += 1
		else:
			read.set_tag(tag, 0)
			ua_counter += 1

		outfile.write(read)

		if (reads_counter % 100000 == 0 and verbose):
			print "##", reads_counter
 
	# Write stats file
	if report:
		handle_stat = open(re.sub(r'\.bam$|\.sam$', '.allelstat', output), 'w')
		handle_stat.write("## " + __file__ + "\n")
		handle_stat.write("## ibam=" + mappedReadsFile + "\n")
		handle_stat.write("## snpFile=" + snpFile + "\n")
		handle_stat.write("## tag=" + tag + "\n")
		handle_stat.write("## output=" + output + "\n")
		handle_stat.write("## verbose=" + str(verbose) + "\n")
		handle_stat.write("## =========================\n")

		handle_stat.write("Total number of snps loaded\t" + str(len(snps)/2) + "\n")
				
		handle_stat.write("## =========================\n")

		handle_stat.write("Total number of reads\t" + str(reads_counter) + "\t100" + "\n")
		handle_stat.write("Number of reads with at least one SNP\t" + str(SNP_counter) + "\t" + str(round(float(SNP_counter)/int(reads_counter)*100,3)) + "\n")
		handle_stat.write("Number of reads assigned to first genome\t" + str(g1_counter) + "\t" + str(round(float(g1_counter)/int(reads_counter)*100,3)) + "\n")
		handle_stat.write("Number of reads assigned to second genome\t" + str(g2_counter) + "\t" + str(round(float(g2_counter)/int(reads_counter)*100,3)) + "\n")
		handle_stat.write("Number of conflicting reads\t" + str(cf_counter) + "\t" + str(round(float(cf_counter)/int(reads_counter)*100,3)) + "\n")
		handle_stat.write("Number of unassigned reads\t" + str(ua_counter) + "\t" + str(round(float(ua_counter)/int(reads_counter)*100,3)) + "\n")
		handle_stat.close()
				
	# Close handler
	infile.close()
	outfile.close()
