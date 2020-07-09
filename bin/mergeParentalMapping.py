#!/usr/bin/env python
import argparse
import os
import pysam
import re
import sys
import time

def argsParse():
    """
    Parsing input & outputs BAM files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paternal", help="Paternal genome alignment file (.bam)")
    parser.add_argument("-m", "--maternal", help="Maternal genome alignment file (.bam)")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("-t", "--tag", help="Allele specific tag", default="XX")
    parser.add_argument("-c", "--comp", help="How to compare alignments (mapq, AS, NM)", default="AS")
    parser.add_argument("--debug",  help="Activate debug mode", action="store_true")

    args = parser.parse_args()
    return (args)


"""
Remove everything after "/" or " " in read's name
"""
def get_read_name(read):
    name = read.query_name
    return re.split('/| ', name)[0]


"""
Get min score on a list of alignment
"""
def getMinScore(reads, score="AS"):
    minScore=100000

    for r in reads:
        if score == "mapq":
            s = r.mapping_quality
        else:
            s = r.get_tag(score)
        if s < minScore:
            minScore = s
        return minScore

"""
Add score to a list of alignment
"""
def setTags(x, key, val):
    for r in x:
        r.set_tag(key, val)
    return(x)

"""
Split a list of aligment into R1/R2 reads
"""
def splitPE(x):
    R1=list()
    R2=list()
    for r in x:
        if r.is_read1:
            R1.append(r)
        elif r.is_read2:
            R2.append(r)
    return R1,R2

"""
Fix a status G1/G2/UA/CF from alignment score
"""
def getAllelicStatus(set1, set2, score="AS", debug=False):
    flag=None
    scoreG1 = getMinScore(set1, score=score)
    scoreG2 = getMinScore(set2, score=score)

    if scoreG1 is not None and scoreG2 is None:
        flag="G1"
    elif scoreG1 is None and scoreG2 is not None:
        flag="G2"
    elif scoreG1 is None and scoreG2 is None:
        flag="UN"
    elif scoreG1 > scoreG2:
        flag="G1"
    elif scoreG1 < scoreG2:
        flag="G2"
    else:
        ## Get all mapping position
        m1pos = list()
        m2pos = list()
        for r in set1:
            m1pos.append(str(r.reference_name) + '_' + str(r.reference_start))
        for r in set2:
            m2pos.append(str(r.reference_name) + '_' + str(r.reference_start))
        if (m1pos == m2pos):
            flag="UA"
        else:
            flag="CF"

    if debug:
        print("[debug] scoreG1=", scoreG1, "- scoreG2=", scoreG2, "- flag=",flag)

    return(flag)

"""
Compare sets of alignments
"""
def compareAlignments(set1, set2, filterNonPrimary=True, asTag="XX", score="AS", debug=False):

    ## Remove unmapped reads
    g1 = [x for x in set1 if not x.is_unmapped]
    g2 = [x for x in set2 if not x.is_unmapped]

    if len(g1) == 0 and len(g2) == 0:
        return(setTags(set1, asTag, "UN"))
    elif len(g1) > 0 and len(g2) == 0:
        return(setTags(set1, asTag, "G1"))
    elif len(g1) == 0 and len(g2) > 0:
        return(setTags(set2, asTag, "G2"))

    ## Remove secondary alignments
    if filterNonPrimary:
        g1 = [x for x in g1 if not x.is_secondary or not x.is_supplementary]
        g2 = [x for x in g2 if not x.is_secondary or not x.is_supplementary]

    flag=getAllelicStatus(g1, g2, score=score, debug=debug)
    if flag == "G2":
        return(setTags(set2, asTag, flag))
    else:
        return(setTags(set1, asTag, flag))


def updateStats(counter_alignments, flag):

    if flag == 'G1':
        counter_alignments['G1'] += 1
    elif flag == 'G2':
        counter_alignments['G2'] += 1
    elif flag == 'UA':
        counter_alignments['UA'] += 1
    elif flag == 'CF':
        counter_alignments['CF'] += 1
    elif flag == 'UN':
        counter_alignments['unmapped'] += 1                                                                                                                                                  

    return(counter_alignments)

"""
Compare two BAM files
"""
def compareBams(bam1, bam2, obam, asTag="XX", score="AS", debug=False):
    
    # Set up counters for report                                                                                                                                                                            
    counter_alignments = {}
    counter_alignments["unmapped"] = 0
    counter_alignments["UA"] = 0
    counter_alignments["G1"] = 0
    counter_alignments["G2"] = 0
    counter_alignments["CF"] = 0
    counter_alignments["total"] = 0

    startInit = None
    isPE = False
    g1Aln=[]
    g2Aln=[]

    with pysam.Samfile(bam1, "rb") as hg1,  pysam.Samfile(bam1, "rb") as hg1Mm, pysam.Samfile(bam2, "rb") as hg2, pysam.Samfile(bam2, "rb") as hg2Mm:
        if hg1.header['SQ'] != hg2.header['SQ']:
            print("SQ header from bam files look different. Stop")
            sys.exit(1)

        out = pysam.AlignmentFile(obam, "wb", template=hg1)
        for g1, g1Mm, g2, g2Mm in zip(hg1.fetch(until_eof=True), hg1Mm.fetch(until_eof=True), hg2.fetch(until_eof=True), hg2Mm.fetch(until_eof=True)):

            ## First initialization - iterate Mm + 1
            if startInit is None:
                g1Mm = next(hg1Mm)
                g2Mm = next(hg2Mm)
                startInit=1
                ## Check if data are paired
                if g1.is_paired and g2.is_paired:
                    isPE=True
                    if debug:
                        print ("[debug] Data are paired")

            if get_read_name(g1) == get_read_name(g2):
                counter_alignments['total'] += 1
                g1Aln.append(g1)
                g2Aln.append(g2)

                ## Deal with multi-Hits and/or PE data
                ## push in g1Aln/g2Aln until the read name changes
                while True:
                    if get_read_name(g1) == get_read_name(g1Mm):
                        g1Aln.append(g1Mm)
                        try:
                            g1 = next(hg1)
                            g1Mm = next(hg1Mm)
                        except StopIteration:
                            break
                    else:
                        break

                while True:
                    if get_read_name(g2) == get_read_name(g2Mm):
                        g2Aln.append(g2Mm)
                        try:
                            g2 = next(hg2)
                            g2Mm = next(hg2Mm)
                        except StopIteration:
                            break
                    else:
                        break

                ## g1/g2 are now on the last Multi reads
                ## g1Mm/g2Mm are on the next read name

                ## Compare two sets of alignments
                if debug:
                    print("[debug] ------------------------")
                    for r in g1Aln: 
                        print("[debug] G1:", r.query_name, r.reference_name, r.reference_start)
                    for r in g2Aln:
                        print("[debug] G2:", r.query_name, r.reference_name, r.reference_start)

                if isPE:
                    counter_alignments['total'] += 1
                    g1PE = splitPE(g1Aln)
                    g2PE = splitPE(g2Aln)

                    bestAlignR1 = compareAlignments(g1PE[0], g2PE[0], asTag=asTag, score=score, debug=debug)
                    flagR1 = bestAlignR1[0].get_tag(asTag)
                    counter_alignments=updateStats(counter_alignments, flagR1)
                    bestAlignR2 = compareAlignments(g1PE[1], g2PE[1], asTag=asTag, score=score, debug=debug) 
                    flagR2 = bestAlignR2[0].get_tag(asTag)
                    counter_alignments=updateStats(counter_alignments, flagR2)

                    if flagR1 == "UN":
                        bestAlignR1.clear()
                    if flagR2 == "UN":
                        bestAlignR2.clear()

                    flag = flagR1 +"/"+ flagR2
                    finalAln = bestAlignR1 + bestAlignR2
                else:
                    bestAlign = compareAlignments(g1Aln, g2Aln, asTag=asTag, score=score, debug=debug)
                    flag = bestAlign[0].get_tag(asTag)
                    counter_alignments = updateStats(counter_alignments, flag)
                    finalAln = bestAlign
                    if flag == "UN":
                        finalAln.clear()

                if finalAln is not None:
                    for r in finalAln:
                        out.write(r)
                    if debug:
                        print("[debug] FLAG:", flag)
                        print("\n")

                ## Drop list elements
                g1Aln.clear()
                g2Aln.clear()
            
            else:
                print("Check that BAM files have the same read names and are sorted in the same way [", get_read_name(g1), "!=", get_read_name(g2),"]")
                sys.exit(1)       

    ## Last occurence of g1/g2
    if get_read_name(g1) != get_read_name(g1Mm) and get_read_name(g2) != get_read_name(g2Mm):
        g1Aln.append(g1Mm)
        g2Aln.append(g2Mm)

        if isPE:
            g1PE = splitPE(g1Aln)
            g2PE = splitPE(g2Aln)
            bestAlignR1 = compareAlignments(g1PE[0], g2PE[0], asTag=asTag, score=score, debug=debug)
            flagR1 = bestAlignR1[0].get_tag(asTag)
            counter_alignments=updateStats(counter_alignments, flagR1)
            bestAlignR2 = compareAlignments(g1PE[1], g2PE[1], asTag=asTag, score=score, debug=debug)
            flagR2 = bestAlignR2[0].get_tag(asTag)
            counter_alignments=updateStats(counter_alignments, flagR2)
            if flagR1 == "UN":
                bestAlignR1.clear()
            if flagR2 == "UN":
                bestAlignR2.clear()
            finalAln = bestAlignR1 + bestAlignR2     
        else:
            bestAlign = compareAlignments(g1Aln, g2Aln, asTag=asTag, score=score, debug=debug) 
            flag = bestAlign[0].get_tag(asTag)
            counter_alignments = updateStats(counter_alignments, flag)
            finalAln = bestAlign
            if flag == "UN":
                finalAln.clear()

        if finalAln is not None:
            for r in finalAln:
                out.write(r)

            if debug:
                print("[debug] FLAG:", finalAln[0].get_tag(asTag))
                print("\n")
  
    hg1.close()
    hg1Mm.close()
    hg2.close()
    hg2Mm.close()

    logName = 'mergeAlignReport.log'
    with open(logName, 'w') as report:
        localtime = time.asctime( time.localtime(time.time()) )

        report.write("------------------------------------------------ \n")
        report.write("| mergeAlign report - " + localtime + " | \n")
        report.write("------------------------------------------------ \n\n")
        report.write("Input parameters\n")
        report.write("================\n")
        report.write("Mapping type:     \tParental\n")
        report.write("Paternal BAM (G1):\t" + bam1 + "\n")
        report.write("Maternal BAM (G2):\t" + bam2 + "\n")
        if not isPE: 
            report.write("Detected read type:\tSingle-end\n")
        if isPE: 
            report.write("Detected reads type:\tPaired-end\n")
        report.write("Output flagged bam:\t" + obam + "\n")
        report.write("\nAllele specific selection report\n")
        report.write("================================\n")
        report.write("Reads specific to G1:                           \t" + str(counter_alignments["G1"]) + \
                     " (" + str(round(float(counter_alignments["G1"])/counter_alignments["total"]*100,3)) + "%)\n")
        report.write("Reads specific to G2:                           \t" + str(counter_alignments["G2"]) + \
                     " (" + str(round(float(counter_alignments["G2"])/counter_alignments["total"]*100,3)) + "%)\n")
        report.write("Reads with no allelic information (UA):         \t" + str(counter_alignments["UA"]) + \
                     " (" + str(round(float(counter_alignments["UA"])/counter_alignments["total"]*100,3)) + "%)\n")
        report.write("Reads with conflicting allelic information (CF):\t" + str(counter_alignments["CF"]) + \
                     " (" + str(round(float(counter_alignments["CF"])/counter_alignments["total"]*100,3)) + "%)\n")
        report.write("Reads unmapped:                                 \t" + str(counter_alignments["unmapped"]) + \
                     " (" + str(round(float(counter_alignments["unmapped"])/counter_alignments["total"]*100,3)) + "%)\n")
        report.write("Total number of reads:                          \t" + str(counter_alignments["total"]) + \
                     " (100%)\n")
        
if __name__ == '__main__':
    args = argsParse()
    compareBams(bam1=args.paternal, bam2=args.maternal, obam=args.output, 
                asTag=args.tag, score=args.comp, debug=args.debug)
