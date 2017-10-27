#!/usr/bin/python

## HiC-Pro
## Copyright (c) 2017 Institut Curie                               
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## The pileup part is tooooooo long. Took 25 hours for 20.5M SNPs
##


"""
Script to count allelic reads for a list of SNPs position. In theory, this script can be replacer by the CheckVariant utils for instance.
Note that the VCF file is loaded in RAM.
"""

import getopt
import sys
import os
import re
import pysam

def usage():
    """Usage function"""
    print "Usage : python countAllelicStatus.py"
    print "-i/--ibam <BAM/SAM file of mapped reads>"
    print "-s/--snp <SNP file information - VCF format>"
    print "[-v/--verbose] <verbose>"
    print "[-h/--help] <Help>"
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:s:vh",
            ["ibam=",
             "snp=",
             "verbose", "help"])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(-1)
    return opts

def get_acgt(samfile, chromosome, start, end):
    # BED files are zero-based as Intervals objects                                                                              
    dcounts = {'A':0, 'C':0, 'G':0, 'T':0, 'tot':0}
    pu = samfile.pileup(chromosome, int(start), int(end))
    for pileupcolumn in samfile.pileup(chromosome, int(start), int(end)):
        if pileupcolumn.pos == int(start):
            for pileupread in pileupcolumn.pileups:
                    #print pileupread.alignment.query_name, pileupread.query_position
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base.upper() in dcounts.keys():
                        dcounts['tot']+=1
                        dcounts[base.upper()]+=1
    return(dcounts)

def is_hetero(basecounts, ref, alt):
    tra = basecounts[ref] + basecounts[alt]

    if tra == 0: return False

    r=float(basecounts[ref])/float(tra)
    ht=float(tra)/basecounts['tot']

    ##print(r, ht)

    if r >= 0.2 and r <= 0.8 and ht > 0.8:
        return True
    else:
        
        return False

mindepth=5

if __name__ == "__main__":

    # Read command line arguments
    opts = get_args()
    inputFile = None
    verbose = False
    
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
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode                                                                                                                                                       
    if verbose:
        print >> sys.stderr, "## " + __file__
        print >> sys.stderr, "## ibam=", mappedReadsFile
        print >> sys.stderr, "## snpFile=", snpFile
        print >> sys.stderr, "## verbose=", verbose, "\n"

    if verbose:
        print >> sys.stderr, "## Opening SAM/BAM file '", mappedReadsFile, "'..."
    samfile = pysam.AlignmentFile(mappedReadsFile, "rb" )

    if verbose:
        print >> sys.stderr, "## Counting allelic information ..."
    vcf_handle = open(snpFile)
    var_counter = 0

    for line in vcf_handle:
        line = line.rstrip()
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            header  = header = line.split('\t')
            header[0]    = header[0][1:]
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
            genotypes  = fields[9].split('\t') if fields[9] else []
            
            ## Counts base at snps position
            basecounts=get_acgt(samfile, chrom, start, start+1)
            
            ## return variants if it is covered enough and is heterozygote
            if basecounts['tot'] >= mindepth and is_hetero(basecounts, ref, alt) :
                print str(chrom) + "\t" + str(start+1) + "\t" + str(ref) + "," + str(alt) + "\t" + str(basecounts[ref]) + "," + str(basecounts[alt]) + "\t" + str(basecounts['tot'])

            ##if basecounts['tot'] >= mindepth and not is_hetero(basecounts, ref, alt) :
            ##    print(chrom, start+1, basecounts['tot'], ref, alt, basecounts, "REJECTED")

        ## verbose 
        if (var_counter % 100000 == 0 and verbose):
            print >> sys.stderr, var_counter

    vcf_handle.close()
    samfile.close()

    
  
    
