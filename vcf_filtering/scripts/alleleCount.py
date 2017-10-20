#!/usr/bin/env python
# -*- coding: utf8 -*-

## Copyright (c) 2016 Institut Curie
## This software is distributed without any guarantee.
## See the LICENCE file for details
                             
## Author(s): Kenzo-Hugo Hillion & Laurène Syx
## Contact(s): kenzo.hillion@curie.fr & laurene.syx@curie.fr
## Python version: 2.7.9
## Script description: Create an allele-specific count table 

scriptVersion = '0.1 - 06-17-2016'

"""
This script is used to create an allele-specific count table 

    INPUTS :  
        - Output file(s) (2 files in case of paired-end sequencing: R1 and R2) from check_variants.py (clinTools) intersected with a BED6 file 
        - VCF file
        - Parent1 and Parent2 names
        
    OUTPUT :
        - Allele-specific count table
    
    USE :
        ./alleleCount.py [opts: --r2 clinToolsFile2 --strand 0:unstranded(default)/1:R1Reverse/2:R1Forward] --r1 clinToolsFile1 -v vcfFile -1 parent1 -2 parent2
"""

###########  Import  ###########

import os
import argparse
import subprocess
import time
import sys
import re
import pandas as pd

###########  Function  ###########

def testNone(argument):
    """
    Test if argument is None or not
    Input(s):
        argument: [int, string, ...] - argument given by the user
    Output:
        [int, string, "", ...]
    """

    if not argument is None:
        variable = argument
    else:
        variable = ""
    return variable

def alleleCount(r1, strand):
    with open(r1, "r") as clinFile:
        snpCount = dict()
        for clinLine in clinFile: 
            if clinLine.startswith('#'):
                pass  
            else:
                baseInfo = clinLine.rstrip().split("\t")
                curChrm = baseInfo[0]
                curStart = baseInfo[1]
                curEnd = baseInfo[2]              
                if strand == 1: 
                    curGene = baseInfo[16]
                    curStrand = baseInfo[18]
                    snpID = str(curChrm) + "_" + str(curStart)+ "." + str(curGene) 
                    if curStrand == "+":
                        countA = baseInfo[6]
                        countT = baseInfo[8]
                        countC = baseInfo[10]
                        countG = baseInfo[12]
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}
                    else:
                        countA = baseInfo[5]
                        countT = baseInfo[7]
                        countC = baseInfo[9]
                        countG = baseInfo[11]
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}     
                elif strand == 2:
                    curGene = baseInfo[16]
                    curStrand = baseInfo[18]
                    snpID = str(curChrm) + "_" + str(curStart)+ "." + str(curGene) 
                    if curStrand == "+":
                        countA = baseInfo[5]
                        countT = baseInfo[7]
                        countC = baseInfo[9]
                        countG = baseInfo[11]
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}
                    else:
                        countA = baseInfo[6]
                        countT = baseInfo[8]
                        countC = baseInfo[10]
                        countG = baseInfo[12]
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}     
                else:
                    curGene = baseInfo[12]
                    curStrand = baseInfo[14]
                    snpID = str(curChrm) + "_" + str(curStart)+ "." + str(curGene) 
                    countA = baseInfo[5]
                    countT = baseInfo[6]
                    countC = baseInfo[7]
                    countG = baseInfo[8]
                    snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}
        return snpCount                

def alleleCountStrandPaired(r1, r2, strand, outname, output_dir):
    file1 = pd.read_csv(r1,sep="\t")
    file2 = pd.read_csv(r2,sep="\t")
    mergeFile = pd.merge(file1, file2, on=["#chromosome", "start", "end", "gene_name", "score", "strand"], how="outer")
    mergeFile = mergeFile.fillna(0)
    mergeFile.to_csv(str(output_dir) + str(outname) + "_mergeCovGene.data",sep="\t",index=False)
    snpCount = dict()
    with open(str(output_dir) + str(outname) + "_mergeCovGene.data", "r") as mergeCovFile:
        for mergeLine in mergeCovFile: 
            if mergeLine.startswith('#'):
                pass  
            else:            
                baseInfo = mergeLine.rstrip().split("\t")
                curChrm = baseInfo[0]
                curStart = int(float(baseInfo[1]))
                curEnd = int(float(baseInfo[2]))
                curGene = baseInfo[16]
                curStrand = baseInfo[18]
                snpID = str(curChrm) + "_" + str(curStart)+ "." + str(curGene)
                if strand == 1:
                    if curStrand == "+":
                        countA = int(float(baseInfo[6])) + int(float(baseInfo[21]))
                        countT = int(float(baseInfo[8])) + int(float(baseInfo[23]))
                        countC = int(float(baseInfo[10])) + int(float(baseInfo[25]))
                        countG = int(float(baseInfo[12])) + int(float(baseInfo[27]))
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}
                    else:
                        countA = int(float(baseInfo[5])) + int(float(baseInfo[22]))
                        countT = int(float(baseInfo[7])) + int(float(baseInfo[24]))
                        countC = int(float(baseInfo[9])) + int(float(baseInfo[26]))
                        countG = int(float(baseInfo[11])) + int(float(baseInfo[28]))
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}     
                elif strand == 2:
                    if curStrand == "+":
                        countA = int(float(baseInfo[5])) + int(float(baseInfo[22]))
                        countT = int(float(baseInfo[7])) + int(float(baseInfo[24]))
                        countC = int(float(baseInfo[9])) + int(float(baseInfo[26]))
                        countG = int(float(baseInfo[11])) + int(float(baseInfo[28]))
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}  
                    else:
                        countA = int(float(baseInfo[6])) + int(float(baseInfo[21]))
                        countT = int(float(baseInfo[8])) + int(float(baseInfo[23]))
                        countC = int(float(baseInfo[10])) + int(float(baseInfo[25]))
                        countG = int(float(baseInfo[12])) + int(float(baseInfo[27]))
                        snpCount[snpID] = {"Chr":curChrm,"Gene":curGene,"Strand":curStrand,"Start":curStart ,"End":curEnd ,"A":countA, "T":countT, "G":countG, "C":countC}   
        return snpCount
            
def getGenotypes(vcf):
    with open(vcf, "r") as vcfFile:
        vcfDict = dict()
        for vcfLine in vcfFile:
            if vcfLine.startswith('#'):
                pass
            else:
                snpInfo = vcfLine.rstrip().split('\t')
                snpChr = str(snpInfo[0])
                snpStart = int(snpInfo[1])-1
                snpEnd = snpInfo[1]
                snpParent1 = snpInfo[3]
                snpParent2 = snpInfo[4]
                snpID = str(snpChr) + "_" + str(snpStart)
                vcfDict[snpID] = {"parent1":snpParent1, "parent2":snpParent2}
        return vcfDict

## Main
if __name__ == "__main__":

    ## Arguments
    parser = argparse.ArgumentParser(prog = "alleleCount.py", description = 'This script is used to create an allele-specific count table',\
                                     epilog = 'Version '+ scriptVersion)

    parser.add_argument("-f", "--r1", type=str, help="R1_file generated by check_variants.py") 
    parser.add_argument("-r", "--r2", type=str, help="R2_file generated by check_variants.py") 
    parser.add_argument("-v", "--vcf", type=str, action="store", help="VCF file")
    parser.add_argument("-1", "--parent1", type=str, action="store",help="Parent 1")
    parser.add_argument("-2", "--parent2", type=str, action="store",help="Parent 2")
    parser.add_argument("-s", "--strand", type = int, help="0 for unstranded seq (default) \
                        / 1 for R1 from a reverse strand / 2 for R1 from a forward strand")
    parser.add_argument('-n', '--name', type = str, \
                        help = 'prefix name for the output file (default : SNPCountTable')
    parser.add_argument('-o', '--output_dir', type = str, \
                        help = 'output directory (default : current directory)')

    args = parser.parse_args()

    ## Store arguments
    r1 = args.r1
    vcf = args.vcf
    parent1 = args.parent1
    parent2 = args.parent2
    
    if (args.r2):
        r2 = args.r2
    else:
        r2 = None

    if (args.name):
        outname = args.name
    else:
        outname = "SNPCountTable"

    if ((args.strand == 0) |(args.strand == 1) | (args.strand == 2)):
        strand = args.strand
    else:
        strand = 0

    if (args.output_dir):
        output_dir = args.output_dir
    else:
        output_dir = "./"
   
    ##########################################################

    ## Allele-specific count (per SNP)
    countTable = open(str(output_dir) + str(outname) + ".data", "w")
    countTable.write("\t".join(map(str,["chromosome", "start", "end", "gene","strand","total", str(parent1)+"_base", str(parent1)+"_count", str(parent2)+"_base", str(parent2)+"_count", "otherCount", "otherBase"]))+"\n")

    if r2 is None:    #If unstranded sequencing or if only stranded single-end sequencing 
        snpCount = alleleCount(r1, strand)
        vcfDict = getGenotypes(vcf)
    else: #If stranded paired-end sequencing
        snpCount = alleleCountStrandPaired(r1, r2, strand,outname, output_dir)
        vcfDict = getGenotypes(vcf)

    for ID, value in zip(snpCount.keys(), snpCount.values()):
        snpID = ID.split(".")[0]
        if snpID in vcfDict:
            snpParent1 = vcfDict[snpID]["parent1"]
            snpParent2 = vcfDict[snpID]["parent2"]
            total = int(value["A"])+int(value["T"])+int(value["C"])+int(value["G"])
            parentCoverage = int(value[snpParent1])+int(value[snpParent2])
            other = int(total)-int(parentCoverage)
            otherBase =  list(set(["A","T","C","G"])-set([snpParent1,snpParent2]))  
            otherCount = ",".join(map(str,[otherBase[0]+":"+str(value[otherBase[0]]), otherBase[1]+":"+str(value[otherBase[1]])]))
            countTable.write("\t".join([str(value["Chr"]), str(value["Start"]),str(value["End"]),str(value["Gene"]), str(value["Strand"]), \
                         str(parentCoverage), str(snpParent1), str(value[snpParent1]), str(snpParent2), str(value[snpParent2]), str(other), str(otherCount)])+"\n")
    countTable.close()
