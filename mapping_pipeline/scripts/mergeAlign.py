#! /usr/bin/env python


## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details
                             
## Author(s): Kenzo-Hugo Hillion & Laurene Syx
## Contact(s): kenzo.hillion@curie.fr & laurene.syx@curie.fr
## Python version: 2.7.9
## Script description: Merging alignments (BAM) from mapping to either parental genomes
## (single or paired-end) or diploid genome (single-end only).

scriptVersion = '1.0 - 09-16-2016'

"""
Script to select alignments based on alignment score ("AS:i") or number of mistmatches ("NM:i")
from either mapping to both parental genomes or a diploid genome.
    
    INPUT(s) :  
        BAM files of mapped alignments sorted by name
    OUTPUT(s) :
        BAM file of selected alignments with a new flag added :
        -XL:i:0 -> alignment is not assigned to any parental genome
        -XL:i:1 -> alignment is assigned to one parental genome
        -XL:i:2 -> alignment is ambiguous
    
    USE :

"""

###########  Import  ###########

import getopt
import sys 
import os
import re
import pysam
import itertools
import random

###########  Constant(s)  ###########

# TEMPORARY FILES
random.seed()
TMP_ID = str(random.randint(1,800000))
MERGED_BAM = TMP_ID + "_parental_merged"
NSORTED_BAM = TMP_ID + "_nsort"

# OUTPUT PARAMETERS/NAMES
FINAL_BAM = "final"
AMBIGUOUS_BAM = "ambiguous"
UNMAPPED_BAM = "unmapped"

# GLOBAL VARIABLES
global AS_TAG,PAR_TAG,WRITE_AMBI,WRITE_UNMAPPED
AS_TAG = "XL"
PAR_TAG = "XP"
WRITE_AMBI = False
WRITE_UNMAPPED = False   

###########  Function(s)  ###########

def usage():
    """Usage function"""
    print "==== mergeAlign.py. version " + scriptVersion + " ===="
    print "Usage : python mergeAlign.py"
    print "For analysis from diploid mapping:"
    print "  -d/--diploid <BAM file from diploid mapping>"
    print "For analysis from parental mapping:"
    print "  -p/--paternal <BAM file from paternal mapping>"
    print "  -m/--maternal <BAM file from maternal mapping>"
    print "Common parameters:"
    print "  -c/--comparison <Method used for the comparison. (1) for Alignment Score and (2) for Number of Mismatches>"
    print "  -s/--sequencing_type <(1) for Single-end and (2) for Paired-end>"
    print "  [-n/--output_name] <Prefix name for output BAM files. Default is \"selected\">"
    print "  [-o/--output_dir] <Output directory. Default is current directory>"
    print "  [-a/--ambiguous] <Write ambiguous alignments in output_dir/prefix_ambiguous.bam>"
    print "  [-u/--unmapped] <Write unmapped alignments in output_dir/prefix_unmapped.bam>"
    print "  [-h/--help] <Help>"
    return

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "d:p:m:c:s:n:o:auh",
            ["diploid=",
             "paternal=",
             "maternal=",
             "comparison=",
             "sequencing_type",
             "output_name",
             "output_dir",
             "ambiguous", "unmapped", "help"])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(-1)
    return opts

def write_single_alignment(alignment,all_outbams):
    '''
    Write the alignment in the correct outbam file

    alignment: alignment [PYSAM object]
    all_outbams: all the different output BAM files [dict]
    '''
    if alignment is None: # Case for paired-end where only one mate is mapped
        return
    elif (alignment.flag >= 512): # Failed quality checks or duplicate, consider unmapped
        if WRITE_UNMAPPED:
            all_outbams["unmapped"].write(alignment)
    elif ((alignment.flag >= 4) & (str(bin(alignment.flag))[-3] == '1')): # unmapped alignment -> 0b0100
        if WRITE_UNMAPPED:
            all_outbams["unmapped"].write(alignment)
    elif (alignment.get_tag(AS_TAG) != 2): # mapped alignment
        all_outbams["final"].write(alignment)
    elif WRITE_AMBI: # ambiguous alignment
        all_outbams["ambi"].write(alignment)
    return

def calc_scores(list_alignments,comparison):
    '''
    This function calculates scores of different alignments from a list
    The smaller the better:
      AS : Alignment Score (-c/--comparison 1)
       -> take absolute value
      NM : Number of mismatch (-c/--comparison 2)

    list_alignments: alignments to calculate score from [list]
    comparison: chosen parameter to compare (AS or NM) [int]

    RETURN scores associated to alignments (same index) [list]
    '''
    scores = []
    for alignment in list_alignments:
            if (alignment.flag >= 512):
                scores.append(1000)
            elif ((alignment.flag >= 4) & (str(bin(alignment.flag))[-3] == '1')):
                scores.append(1000)
            elif (comparison == 1):
                scores.append(abs(alignment.get_tag("AS")))
            elif (comparison == 2):
                scores.append(alignment.get_tag("NM"))
    return scores


def comp_two_single_alignments(list_alignments,scores):
    '''
    This function compares two alignments on scores where smaller is better
    and returns the best alignment with a TAG (AS_TAG):
        0: Equally mapped on both parents
        1: Preferentially mapped on one parent
        2: Ambiguous

    list_alignments: 2 alignments to compare [list]
    scores: 2 scores associated to alignments (same index) [list]
    '''
    alignment1 = list_alignments[0]
    alignment2 = list_alignments[1]
 
    # If both alignments unmapped
    if ((scores[0] == 1000) & (scores[1] == 1000)):
        alignment1.set_tag(AS_TAG,0)
        return alignment1

    # Reads mapped at the same position on the genome
    if (scores[0] < scores[1]): # Best score from parent 1
        alignment1.set_tag(AS_TAG,1)
        return alignment1
    elif (scores[0] > scores[1]): # Best score from parent 2
        alignment2.set_tag(AS_TAG,1)
        return alignment2
    else: # Same score
        if ((alignment1.pos == alignment2.pos) & (alignment1.rname == alignment2.rname)):
            # No allelic information on this read
            alignment1.set_tag(AS_TAG,0)
            return alignment1
        else:
            # Ambiguous case
            alignment1.set_tag(AS_TAG,2)
            return alignment1 

def comp_three_single_alignments(list_alignments,scores):
    '''
    This function compares three alignments on scores where smaller is better
    and returns the best alignment with a TAG (AS_TAG):
        0: Equally mapped on both parents
        1: Preferentially mapped on one parent
        2: Ambiguous

    list_alignments: 3 alignments to compare [list]
    scores: 3 scores associated to alignments (same index) [list]
    '''
    
    if (max(scores) == min(scores)): # Ambiguous case : at least 2 positions with same score
        alignment = list_alignments[0]
        alignment.set_tag(AS_TAG,2)
        return alignment
    else:
        # Get first minimum score
        alignment1 = list_alignments[scores.index(min(scores))]
        score1 = scores[scores.index(min(scores))]
        # Remove alignment and score from lists
        del list_alignments[scores.index(min(scores))]
        del scores[scores.index(min(scores))]
        alignment2 = list_alignments[scores.index(min(scores))]
        score2 = scores[scores.index(min(scores))]
        return comp_two_single_alignments([alignment1,alignment2],[score1,score2])

def comp_single_alignments(list_alignments,comparison):
    '''
    Function to select and return the best position for a alignment either:
      mapped on both parental genomes
      mapped on a diploid genome
      
    list_alignments: different lines for the alignment [list]
    comparison: parameter chosen for the comparison [int]
    '''
    scores = calc_scores(list_alignments,comparison)
    if (len(list_alignments) > 2):
    # Three alignments have been reported
        alignment = comp_three_single_alignments(list_alignments,scores)
    elif (len(list_alignments) > 1):
    # If only two alignments were reported
        alignment = comp_two_single_alignments(list_alignments,scores)
    else:
        alignment = list_alignments[0]
        if ((alignment.flag >= 4) & (str(bin(alignment.flag))[-3] == '1')):
            alignment.set_tag(AS_TAG,0)
        else:
            alignment.set_tag(AS_TAG,1)
    return alignment

def init_pair(pairs,tag):
    '''
    function to initialize a dict entry

    pairs: all alignments pairs for a pair of read [dict]
    tag: genotype of origin of alignment stored in PAR_TAG [string]
    '''
    pairs[tag] = {}
    pairs[tag]["R1"] = None
    pairs[tag]["R2"] = None
    pairs[tag]["score"] = 0
    pairs[tag]["mapped"] = 0 # how many mapped from the pair ?
    return pairs

def reformat_pairs(list_alignments,scores):
    '''
    Reformat all alignments for a pair of read in a dictionary
    
    list_alignments: different lines for the alignment and its mate [list]
    scores: scores corresponding to each line [list]

    WARNING: strategy only working for parental mapping so far...

    RETURN: dictionary of position [dict]:
      [PAR_TAG] --->  ["R1"] ---> First in pair alignment (Def: None)
                      ["R2"] ---> Second in pair alignment (Def: None)
                      ["score"] ---> Addition of both scores for comparison
                      ["mapped"] ---> how many mapped from the pair
    '''
    pairs = {}
    for i in range(0,len(list_alignments)):
        # Unstack from both list
        alignment = list_alignments.pop()
        score = scores.pop()
        
        # First or second in pair ? 
        if ((alignment.flag >= 128) & (str(bin(alignment.flag))[-8] == '1')):
            # This is the second in pair
            if not (alignment.get_tag(PAR_TAG) in pairs.keys()):
                pairs = init_pair(pairs,alignment.get_tag(PAR_TAG))
            pairs[alignment.get_tag(PAR_TAG)]["R2"] = alignment
        elif ((alignment.flag >= 64) & (str(bin(alignment.flag))[-7] == '1')):
            # This is the first in pair
            if not (alignment.get_tag(PAR_TAG) in pairs.keys()):
                pairs = init_pair(pairs,alignment.get_tag(PAR_TAG))
            pairs[alignment.get_tag(PAR_TAG)]["R1"] = alignment
        else:
            # Flag < 64 read not paired
            sys.stderr.write("WARNING: " + alignment.qname + " is not paired. Skipped")
            continue

        # Increase score and mapped counter
        pairs[alignment.get_tag(PAR_TAG)]["score"] += score
        if not (str(bin(alignment.flag))[-3] == '1'):
            pairs[alignment.get_tag(PAR_TAG)]["mapped"] += 1
    return pairs

def tag_selected_pair(pair,tag):
    '''
    Function that tags a selected pair of alignment
      0 if not mapped
      1 if mapped, meaning that it has allelic information

    pair: list alignments for the pair [list] of size 2
    tag: tag to assign [int]
    '''
    for alignment in pair:
        if alignment is not None:
            if ((alignment.flag >= 4) & (str(bin(alignment.flag))[-3] == '1')):
                # Read not mapped
                alignment.set_tag(AS_TAG,0)
            else:
                alignment.set_tag(AS_TAG,tag)
    return pair

def comp_two_paired_alignments(dict_pair):
    '''
    Function to select the best pair of position between two pairs.

    dict_pair: two pairs for each parent [dict]

    RETURN: Best pair [list] os size 2
    '''
    genos = dict_pair.keys()
    if (dict_pair[genos[0]]["mapped"] > dict_pair[genos[1]]["mapped"]):
        # more reads aligned in genos[0] VS genos[1]
        # Keep alignments from genos[0]
        pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
        tag_selected_pair(pair,1)
    elif (dict_pair[genos[0]]["mapped"] < dict_pair[genos[1]]["mapped"]):
        # less reads aligned in genos[0] VS genos[1]
        # Keep alignments from genos[1]
        pair = [dict_pair[genos[1]]["R1"],dict_pair[genos[1]]["R2"]]
        tag_selected_pair(pair,1)
    else:
        # Same number of reads mapped in both genotypes (either 1 or 2)
        # Need to compare both pairs on scores
        if (dict_pair[genos[0]]["score"] < dict_pair[genos[1]]["score"]):
            # Score better for genos[0]
            pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
            tag_selected_pair(pair,1)
        elif (dict_pair[genos[0]]["score"] > dict_pair[genos[1]]["score"]):
            # Score better for genos[1]
            pair = [dict_pair[genos[1]]["R1"],dict_pair[genos[1]]["R2"]]
            tag_selected_pair(pair,1)
        else:
            # Both pairs have the same score
            if (dict_pair[genos[0]]["mapped"] == 2):
                # Both reads are mapped
                if ((dict_pair[genos[0]]["R1"].pos == dict_pair[genos[1]]["R1"].pos) & \
                (dict_pair[genos[0]]["R1"].rname == dict_pair[genos[1]]["R1"].rname) & \
                (dict_pair[genos[0]]["R1"].pnext == dict_pair[genos[1]]["R1"].pnext)):
                    # Same positions for both pairs
                    pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                    tag_selected_pair(pair,0)
                else:
                    # Ambiguous case, different positions with same score
                    pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                    tag_selected_pair(pair,2)
            else:
                # Only one read has been mapped in both case
                if ((dict_pair[genos[0]]["R1"] is not None) & (dict_pair[genos[1]]["R1"] is not None)):
                    # R1 mapped for both
                    if ((dict_pair[genos[0]]["R1"].pos == dict_pair[genos[1]]["R1"].pos) & \
                    (dict_pair[genos[0]]["R1"].rname == dict_pair[genos[1]]["R1"].rname)):
                        # Same positions for both reads
                        pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                        tag_selected_pair(pair,0)
                    else:
                        # Ambiguous case, different positions with same score
                        pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                        tag_selected_pair(pair,2)
                elif ((dict_pair[genos[0]]["R2"] is not None) & (dict_pair[genos[1]]["R2"] is not None)):
                    # R2 mapped for both
                    if ((dict_pair[genos[0]]["R2"].pos == dict_pair[genos[1]]["R2"].pos) & \
                    (dict_pair[genos[0]]["R2"].rname == dict_pair[genos[1]]["R2"].rname)):
                        # Same positions for both reads
                        pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                        tag_selected_pair(pair,0)
                    else:
                        # Ambiguous case, different positions with same score
                        pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                        tag_selected_pair(pair,2)
                else:
                    # Ambiguous cases
                    pair = [dict_pair[genos[0]]["R1"],dict_pair[genos[0]]["R2"]]
                    tag_selected_pair(pair,2)
    return pair

def comp_three_paired_alignments(dict_pair):
    '''
    
    '''
    sys.stderr.write("WARNING: more than 2 pair of alignment detected. Skipped")
    return [None,None]

def comp_paired_alignments(list_alignments,comparison):
    ''' 
    Function to select the best position for paired-end alignments either:
      mapped on both parental genomes
      [NOT IMPLEMENTED] : mapped on a diploid genome
      
    list_alignments: different lines for the alignment and its mate [list]
    comparison: parameter chosen for the comparison [int]
    
    WARNING: strategy only working for parental mapping so far...

    RETURN: best pair [list] of size 2
    '''
    scores = calc_scores(list_alignments,comparison)
    # Reformat alignments in a dictionnary:
    dict_pair = reformat_pairs(list_alignments,scores)
    if (len(dict_pair) > 2):
    # Three alignments have been reported
        pair = comp_three_paired_alignments(dict_pair)
    elif (len(dict_pair) > 1):
    # Alignments from both parents
        pair = comp_two_paired_alignments(dict_pair)
    elif (len(dict_pair) > 0):
    # Alignments only from one parent
        pair = [dict_pair[dict_pair.keys()[0]]["R1"],dict_pair[dict_pair.keys()[0]]["R2"]]
        tag_selected_pair(pair,1)
    else:
    # No pair recorded
        pair = [None,None]
    return pair

def create_header(header):
    '''
    Function to create a BAM header from a diploid BAM file
      to remove the repetition of each chromosome from diploid genome
      ex: chrX_CAST and chrX_129 -> chrX

    header: header to be modified [dict]
    
    RETURN : dictionnary of correspond ID for each chromosome [dict]
    '''
    chrom_ids = {} # IDs of chromosomes to be return
    new_sq_header = [] # new header to replace old one
    cpt = 0
    for ref in header['SQ']:
        chrom = ref['SN'].split("_",1)[0]
        # Skip if chromosome already exists
        if chrom in chrom_ids.keys():
            continue
        # Build new chromosome information
        tmp_dic = {}
        tmp_dic['LN'] = ref['LN']
        tmp_dic['SN'] = chrom
        # Add new chromosome information to the header
        new_sq_header.append(tmp_dic)
        # Save the position of the chromosome in the SQ header part
        chrom_ids[chrom] = cpt
        cpt += 1
    # Replace the SQ part of the header by the new one
    header['SQ'] = new_sq_header
    return chrom_ids

def rephase_alignment(alignment,chrom_ids):
    '''
    Change reference_id of the alignment to the new corresponding one
    from a diploid mapping BAM file
      Also add the information in PAR_TAG
      from which genotype the alignment comes from

    alignment: alignment [PYSAM object]
    chrom_ids: corresponding IDs for each chromosome [dict]
    '''
    # Check if read unmapped, just skip
    if ((alignment.flag >= 4) & (str(bin(alignment.flag))[-3] == '1')):
        return
    # Get the name of the chromosome
    chrom_id = chrom_ids[alignment.reference_name.split("_",1)[0]]
    # Save origin of the chromosome in PAR_TAG tag
    alignment.set_tag(PAR_TAG,alignment.reference_name.split("_",1)[1])
    # Replace with the new ID corresponding to output header
    alignment.reference_id = chrom_id
    return


##########  Main ##########

if __name__ == "__main__":

    ## Read command line arguments
    opts = get_args()

    ## Default parameters and initialization
    diploid = None
    paternal = None
    maternal = None
    comparison = 1 # Default is 1 for alignment score
    seq_type = 1 # Default is 1 for single end
    outname = "selected_" # Default name is "selected"
    outdir = "./" # Default name is "selected"
    mapping_diploid = False
    mapping_parental = False

    ## Assign command line arguments
    if len(opts) == 0:
        usage()
        sys.exit()

    for opt,arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--diploid"):
            diploid = arg
            mapping_diploid = True
        elif opt in ("-p", "--paternal"):
            paternal = arg
            mapping_parental = True
        elif opt in ("-m", "--maternal"):
            maternal = arg
        elif opt in ("-c", "--comparison"):
            comparison = int(arg)
        elif opt in ("-s", "--sequencing_type"):
            seq_type = int(arg)
        elif opt in ("-o", "--output_dir"):
            outdir = str(arg) + "/"
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        elif opt in ("-n", "--output_name"):
            outname = str(arg) + "_"
        elif opt in ("-a", "--ambiguous"):
            WRITE_AMBI = True
        elif opt in ("-u", "--unmapped"):
            WRITE_UNMAPPED = True
        else:
            assert False, "unhandled option"


    ## Test arguments

    ###########################################################

    ## == Merge (Parental) and sort BAM file ==

    if (mapping_parental):
        # Merge BAM files
        mergedbam = None
        for bam_name in paternal,maternal:
            bam = pysam.Samfile(bam_name, "rb")
            # Output merged file
            if (mergedbam == None):
                mergedbam = pysam.Samfile(outdir + outname + MERGED_BAM + ".bam", "wb", template=bam)
            # Get the name of the bam from which the alignments come from
            origin = bam_name.split(".",1)[0]
            for alignment in bam.fetch(until_eof=True):
                alignment.set_tag(PAR_TAG,origin)
                mergedbam.write(alignment)
            # Closing BAM file
            bam.close()
        # Closing merged BAM
        mergedbam.close()

    # Sort BAM file by name for analysis
    if (mapping_parental):
        pysam.sort("-n",outdir + outname + MERGED_BAM + ".bam",outdir + outname + NSORTED_BAM)
        os.remove(outdir + outname + MERGED_BAM + ".bam")
    elif (mapping_diploid):
        pysam.sort("-n",diploid,outdir + outname + NSORTED_BAM)


    ## == Selecting best position ==

    # Opening name sorted file
    sortedbam = pysam.Samfile(outdir + outname + NSORTED_BAM + ".bam", "rb")
    # Final output file
    if (mapping_parental):
        finalbam = pysam.Samfile(outdir + outname + FINAL_BAM + ".bam", "wb", template=sortedbam)
    elif (mapping_diploid):
        # Declare empty dictionnary to store chrom_IDs of new header
        chrom_ids = {}
        # Get previous header to build new one
        new_bam_header = sortedbam.header
        # Build new header
        chrom_ids = create_header(new_bam_header)
        # Finally open output BAM file with new header
        finalbam = pysam.Samfile(outdir + outname + FINAL_BAM + ".bam", "wb", header=new_bam_header)

    # Dictionnary with all output BAM files
    all_outbams = {}
    all_outbams["final"] = finalbam

    if (WRITE_AMBI):
        print ("Writing ambiguous alignments in: " + outdir + outname + AMBIGUOUS_BAM + ".bam")
        if (mapping_parental):
            ambibam = pysam.Samfile(outdir + outname + AMBIGUOUS_BAM + ".bam", "wb", template=sortedbam)
        elif (mapping_diploid):
            ambibam = pysam.Samfile(outdir + outname + AMBIGUOUS_BAM + ".bam", "wb", header=new_bam_header)
        all_outbams["ambi"] = ambibam
    if (WRITE_UNMAPPED):
        print ("Writing unmapped alignments in: " + outdir + outname + UNMAPPED_BAM + ".bam")
        if (mapping_parental):
            unmapbam = pysam.Samfile(outdir + outname + UNMAPPED_BAM + ".bam", "wb", template=sortedbam)
        elif (mapping_diploid):
            unmapbam = pysam.Samfile(outdir + outname + UNMAPPED_BAM + ".bam", "wb", header=new_bam_header)
        all_outbams["unmapped"] = unmapbam

    # Set up counters
    counter_alignment = 0

    # Variables inititialization
    prev_alignments = []
 
    for alignment in sortedbam.fetch(until_eof=True):
        counter_alignment += 1
        if (mapping_diploid):
            rephase_alignment(alignment,chrom_ids)
        if (len(prev_alignments) == 0):
            prev_alignments.append(alignment)
        else:
            if (alignment.qname == prev_alignments[0].qname):
                prev_alignments.append(alignment)
            else:
                # First we need to process the previous group of alignments
                if (seq_type == 1): # Single-end
                    # Select best position
                    selected_alignment = comp_single_alignments(prev_alignments,comparison)
                    # Write alignment in output file
                    write_single_alignment(selected_alignment,all_outbams)
                elif (seq_type == 2): # Paired-end
                    selected_alignment = comp_paired_alignments(prev_alignments,comparison)
                    for al_pair in selected_alignment:
                        write_single_alignment(al_pair,all_outbams)
                prev_alignments = []
                # Then we process the alignment as it is a first alignment
                prev_alignments.append(alignment)

        if (counter_alignment % 1000000 == 0):
            print ("Reads treated: " + str(counter_alignment))

    # Treat last stack of alignments
    if (len(prev_alignments) > 0):
        # Process last group of alignment
        if (seq_type == 1): # Single-end
            # Select best position
            selected_alignment = comp_single_alignments(prev_alignments,comparison)
            if (mapping_diploid):
                rephase_alignment(selected_alignment,chrom_ids)
            # Write alignment in output file
            write_single_alignment(selected_alignment,all_outbams)
        elif (seq_type == 2): # Paired-end
            selected_alignment = comp_paired_alignments(prev_alignments,comparison)

    print ("Total number of alignments: " + str(counter_alignment))

    # Closing SAM/BAM files
    for key in all_outbams:
        all_outbams[key].close()
    sortedbam.close()
    
    os.remove(outdir + outname + NSORTED_BAM + ".bam")
