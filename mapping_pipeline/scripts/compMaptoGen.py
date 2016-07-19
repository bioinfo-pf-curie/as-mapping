#! /usr/bin/env python


## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details
    
## Author(s): Kenzo-Hugo Hillion
## Contact(s): kenzo.hillion@curie.fr
## Python version: 2.7
## Script description: Script to compare generated BAM (from ART) to aligned BAM.

scriptVersion = '0.1 - 05-24-2016'

"""
Script to compare generated BAM (from ART) to aligned BAM.

    INPUT(s):
        BAM files sorted by names
        Strain names
        Output directory
    OUTPUT(s):
        mapping_status.bam     : BAM file with XM flag (1 : mapped to the correct region,
                                                        2 : mapped to incorrect region
                                                        3 : not mapped)
                                               XP flag : position of origin

        allele_counts.txt    : Counts for each SNPs corresponding to the interval 
                               (if no SNP, correspond to the middle of the interval)
    
    USE :
        compMaptoGen.py -g generated.bam -m mapped.bam -1 geno1 -2 geno2 -o OUTPUT_DIR/

"""

###########  Import  ###########

import sys
import getopt
import os
import pysam
import itertools

###########  Function(s)  ###########

def usage():
    """Usage function"""
    print "==== compMaptoGen.py. version " + scriptVersion + " ===="
    print "To compare BAM file generated by ART VS mapped BAM."
    print " -> Also write an output bam with a flag"
    print "  -XM=1 : mapped desired region"
    print "  -XM=2 : mapped to incorrect region"
    print "  -XM=3 : not mapped"
    print "  -XP : corresponds to the position of origin"
    print ""
    print "Usage : compMaptoGen.py"
    print "-g/--genbam <name sorted BAM file generated by ART (simulation reads tool)>"
    print "-m/--mapbam <name sorted BAM file of mapped reads>"
    print "-1/--strain1 <name of strain 1>"
    print "-2/--strain2 <name of strain 2>"
    print "[-o/--output_dir] <Output directory>"
    print "[-u/--unmapped] <report XM:i:2 and XM:i:3 from mapped and generated BAM in a supplementary BAM file>"
    print "[-h/--help] <Help>"
    return

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
        sys.argv[1:],
        "g:m:o:1:2:uh:",
        ["genbam","mapbam","output_dir","strain1","strain2","unmapped","help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

###########  Constant(s)  ###########

OUTBAM = "mapping_status.bam"
REPBAM = "map_gen_XM2-3.bam"
OUTCOUNTS = "allele_counts.txt"
TAG_MAP = "XM"
TAG_POS = "XP"

###########  Main  ###########

if __name__ == "__main__":
    
    # Script arguments
    opts = get_args()
    if len(opts) == 0:
        usage()
        sys.exit()

    # Default values for arguments    
    genbam = None
    mapbam = None
    rep_xm = False
    strain1 = ""
    strain2 = ""
    output_dir = ""

    # Store arguments
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-g", "--genbam"):
            genbam = arg
        elif opt in ("-m", "--mapbam"):
            mapbam = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-1", "--strain1"):
            strain1 = arg
        elif opt in ("-2", "--strain2"):
            strain2 = arg
        elif opt in ("-u", "--unmapped"):
            rep_xm = True
        else:
            assert False, "unhandled option"
    
    # Test arguments
    if ((strain1 == "") | (strain2 == "")):
        print "INPUT ERROR : strain name(s) not specified in arguments."
        usage()
        sys.exit()
    if ((genbam is None) | (mapbam is None)):
        print "INPUT ERROR : BAM file(s) not specified in arguments."
        usage()
        sys.exit()

    ##################################################################

    # Read SAM/BAM files
    filegenbam = pysam.Samfile(genbam, "rb")
    filemapbam = pysam.Samfile(mapbam, "rb")

    # Output BAM file
    outfile = pysam.AlignmentFile(output_dir + OUTBAM, "wb", template=filemapbam)
    if (rep_xm):
        gen_xm = pysam.AlignmentFile(output_dir + REPBAM, "wb", template=filemapbam)
    
    # Counting correctly/incorrectly mapped reads
    counter_read = 0
    read_length = 0
    dict_counter = {}
    for x,y in itertools.izip(filemapbam.fetch(until_eof=True),filegenbam.fetch(until_eof=True)):
        counter_read += 1                            # Increments read counter
        if (read_length == 0):                        # Calculate read_length
            read_length = len( x.seq )
        tmp_split=x.qname.split("-", 1 )            # Name of the generated reads : chr:begin-end_strain-nread
        pos=int(tmp_split[0].split(":", 1)[1])+read_length    # Get the position of the SNP
        if not(dict_counter.has_key(pos)):
            dict_counter[pos] = {}                    # Create dictionnary and initialize each counter
            dict_counter[pos][strain1] = 0
            dict_counter[pos][strain2] = 0
            dict_counter[pos]["total"] = 0
            dict_counter[pos]["wrong"] = 0
            dict_counter[pos]["unmapped"] = 0
        dict_counter[pos]["total"] += 1                # Increments total read for this interval
        x.set_tag(TAG_POS,str(y.reference_name)+'-'+str(y.pos),value_type="Z")
        if (x.pos==y.pos):                            # Read mapped correctly (consider that read mapping to the exact same region on another chromosome is pretty unlikely)
            tmp_split=x.qname.split("_", 1 )        # Get the name of the strain
            strain=tmp_split[1].split("-",1)[0]    
            dict_counter[pos][strain] += 1            # Increments the counter for the corresponding strain
            x.set_tag(TAG_MAP,1)                        # Set tag XM:i:1 for correctly mapped read
        else:                                        # Read mapped to incorrect region or unmapped
            if (x.flag == 4):                        # Read not mapped
                dict_counter[pos]["unmapped"] += 1    # Increments unmapped counter
                x.set_tag(TAG_MAP,3)                        # Set tag XM:i:3 for read not mapped
            else:
                dict_counter[pos]["wrong"] += 1        # Increments incorrect mapped read counter
                x.set_tag(TAG_MAP,2)                        # Set tag XM:i:2 for incorrectly mapped read
            if (rep_xm):
                gen_xm.write(x)
        outfile.write(x)
        if (counter_read % 1000000 == 0):
            print ("Reads treated : " + str(counter_read))
    print ("Total number of reads treated : " + str(counter_read))

    # Writing output file with counts
    print ("Writing output files in " + str(output_dir))
    COUNTS = open(str(output_dir) + OUTCOUNTS, "w")
    COUNTS.write("SNP\t" + strain1 + "\t" + strain2 + "\tTotal_reads\tUnmapped_reads\tWrongly_mapped\n")
    for read in dict_counter.keys():
        COUNTS.write(str(read) + "\t" + str(dict_counter[read][strain1]) + "\t" \
+ str(dict_counter[read][strain2]) + "\t" + str(dict_counter[read]["total"]) + "\t" \
+ str(dict_counter[read]["unmapped"]) + "\t" + str(dict_counter[read]["wrong"]) + "\n")
    # Closing output file
    COUNTS.close()

    # Closing SAM/BAM files
    filegenbam.close()
    filemapbam.close()
    if (rep_xm):
        gen_xm.close()
    outfile.close()
