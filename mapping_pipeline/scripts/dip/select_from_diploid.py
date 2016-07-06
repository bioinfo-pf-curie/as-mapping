#! /usr/bin/env python


## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details

## Author(s): Kenzo-Hugo Hillion
## Contact(s): kenzo.hillion@curie.fr
## Python version: 2.7
## Script description: Script to select a read based on alignment score ("AS:i") mapped on a diploid genome

scriptVersion = '0.1 - 06-03-2016'

"""
Script to select a read based on alignment score ("AS:i") mapped on a diploid genome
    

    INPUT(s) :    
        BAM file of mapped reads sorted by name (asked bowtie2 to return 2 best reads)
    OUTPUT(s) :
        BAM file of selected reads :
        -XA:i:0 -> read is not assigned to any parental genome
        -XA:i:3 -> read is assigned to any parental genome because of mapping to different region with the same alignment score
        -XA:i:4 -> read is assigned to a position (but mapped to different position with lower score on second genome)   
 
    USE :
        ./select_best.py -b diploid.bam -o output_dir/

"""

###########  Import  ###########

import sys
import getopt
import os
import re
import pysam

########## Function(s) ##########

def usage():
    print "==== select_from_diploid.py. version " + scriptVersion + " ===="
    print "Select read mapped on a diploid genome from the BAM file sorted by names (with 3 best alignments returned)."

    """Usage function"""
    print "Usage : select_from_diploid.py"
    print "-b/--bam <names sorted BAM file of mapped reads (2 best alignments returned)>"
    print "[-n/--name] <prefix name for output bam>"
    print "[-o/--output_dir] <Output directory>"
    print "[-h/--help] <Help>"
    return

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
        sys.argv[1:],
        "h:o:b:n:",
        ["help","output_dir","bam","name"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

def comp_two_reads(read1,read2,tag):
    '''
    This function compares the position and alignment score (AS) of two reads
    and return the read at the position with the best AS 
    or a read tagged with XA:i:3 for ambiguous case
    '''
    if ((read1.pos == read2.pos) & (read1.rname == read2.rname)):
    # Reads mapped on the same chromosome at the same position on both genome
        if (read1.get_tag("AS") >= read2.get_tag("AS")):
            return read1
        else:
            return read2
    else:
    # Reads mapped at 2 different locations
        if (read1.get_tag("AS") == read2.get_tag("AS")):
        # Ambiguous case
            read1.set_tag(tag,3)
            return read1
        elif (read1.get_tag("AS") > read2.get_tag("AS")):
            read1.set_tag(tag,4)
            return read1
        else:  
        # This case should never happen as the reads are supposed to be
        # ranked by Alignment Scores
            read2.set_tag(tag,4)
            return read2


def comp_three_reads(read1,read2,read3,tag):
    '''
    This function compares alignment scores (AS) of three reads
    and returns the read at the position with the best AS
    or a read tagged with XA:i:3 for ambiguous case
    '''
    if (read1.get_tag("AS") > read3.get_tag("AS")):
    # We can only compare the two first reads
        return comp_two_reads(read1,read2,tag)
    else:
    # This case is ambiguous because AS from first read is either equal
    # or greater than the two others.
    # It means here that the three AS are identical
        read1.set_tag(tag,3)
        return read1


def comp_reads(list_reads,tag):
    '''
    Function description
    '''
    if (len(list_reads) > 2):
    # Three reads have been reported
        return comp_three_reads(list_reads[0],list_reads[1],list_reads[2],tag)
    elif (len(list_reads) > 1):
    # If only two reads were reported
        return comp_two_reads(list_reads[0],list_reads[1],tag)
    else:
        return list_reads[0]
    return list_reads[0]

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
    bam=None

    # Store arguments
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-b", "--bam"):
            bam = arg
        elif opt in ("-n", "--name"):
            outname = arg
        else:
            assert False, "unhandled option"

    # Test arguments
    if (bam is None):
        print "INPUT ERROR : BAM file not specified in arguments."
        usage()
        sys.exit()

    ##################################################################

    # Read SAM/BAM files
    bamfile = pysam.Samfile(bam, "rb")

    # Open a file for writing
    bamout = pysam.Samfile(output_dir + outname + ".bam", "wb", template=bamfile)

    # Set up counters and name of the tag
    counter_read = 0
    tag = "XA"

    # Variables initialization
    as_r1 = 0
    as_r2 = 0
    prev_read = []
    
    for read in bamfile.fetch(until_eof=True):
        counter_read += 1        # Increments read counter
        if (len(prev_read) == 0):
            if (read.flag == 4):    # If the read is unmapped
                read.set_tag(tag,0)
                bamout.write(read)
                continue
            prev_read.append(read)  # Save the read for comparison to the next one
        
        else:
            if (read.qname == prev_read[0].qname):
                prev_read.append(read)
            else:
                # First we need to process the previous group of reads
                bamout.write(comp_reads(prev_read,tag))
                prev_read = [] # Reinitiliaze the list

                # Then we process the read as it is a first read
                if (read.flag == 4):    # If the read is unmapped
                    read.set_tag(tag,0)
                    bamout.write(read)
                    continue
                prev_read.append(read)  # Save the read for comparison to the next one

        if (counter_read % 1000000 == 0):
            print ("Reads treated : " + str(counter_read))
    
    if (len(prev_read) > 0):
        bamout.write(comp_reads(prev_read,tag))
    
    print ("Total number of reads : " + str(counter_read))

    # Closing SAM/BAM files
    bamout.close()
    bamfile.close()
