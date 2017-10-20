#! /usr/bin/env python
# -*- coding: utf8 -*-

# This file is part of clinTools
# Copyright (C) 2015 Vivien Deshaies - Institut Curie
# 
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, tempfile, subprocess, sys, re, string
import hgvs.variantmapper
import logging

STDERR_FILE = sys.stderr
TRANSLATION_TAB = string.maketrans('atcgATCG', 'tagcTAGC')

def checkSamtoolsVersion(version):
    sp = subprocess.Popen("samtools", stderr = subprocess.PIPE)
    stderr = str(sp.communicate())
    match = re.findall(r'Version: (.*?) ', stderr)
    if match[0] != version:
        logging.warning("The samtools version in your environnement is {0}, this tool was tested with samtools version {1}\nTo avoid errors you may consider updating your PATH environnement variable or update your samtools version".format(match[0], version))

def parseChrGtf(gtfFile):

    """Parse a gtf file in dict

    :param gtfFile: a gtf file path
    :type gtfFile: str
    :return annotDict: a dictionnary of transcripts annotations
    :rtype: dict
    """

    nmPattern = re.compile(r"transcript_id \"(NM_[0-9]+[_dup.]*[1-9]*)\"")
    genePattern = re.compile(r"gene_id \"(.*?)\"")
    annotDict = dict()
    acceptedNames = ['start_codon', 'stop_codon', 'exon']

    with open(gtfFile, 'r') as in_handle:
        for line in in_handle:
            line = line.split("\t")
            chrom = line[0]
            name = line[2]
            # convert to 0-base index
            start = int(line[3]) -1
            end = int(line[4]) - 1
            strand = line[6]

            if name not in acceptedNames:
                continue

            # get nm id
            nmMatch = nmPattern.search(line[8])
            if nmMatch:
                nm = nmMatch.group(1)
            else:
                continue

            # get gene name
            geneMatch = genePattern.search(line[8])
            if geneMatch:
                gene = geneMatch.group(1)
            else:
                continue

            if chrom not in annotDict.keys():
                annotDict[chrom] = dict()
            
            if nm not in annotDict[chrom].keys():
                annotDict[chrom][nm] = {
                    'start_codon' : list(),
                    'stop_codon' : list(),
                    'exons' : list(),
                    'gene' : gene,
                    'strand' : ""
                }

            rnaDict = annotDict[chrom][nm]

            if name == 'start_codon':
                rnaDict['start_codon'] = [int(start), int(end)]
            elif name == 'stop_codon':
                rnaDict['stop_codon'] = [int(start), int(end)]
            elif name == 'exon':
                rnaDict['exons'].append([int(start), int(end)])
            
            if rnaDict['strand'] == "":
                rnaDict['strand'] = strand
            elif rnaDict['strand'] != strand:
                raise Exception("All features must be on the same strand. nm :" + nm + " line : " + str(line)) 

    for chromDict in annotDict.values():
        for rnaDict in chromDict.values():
            # sort exons by start codon
            if rnaDict['strand'] == "+":
                rnaDict['exons'] = sorted(rnaDict['exons'], key=lambda exon : exon[0]) 

            elif rnaDict['strand'] == "-":
                # revert all start stop list to have it in the proper order
                rnaDict['start_codon'].reverse()
                rnaDict['stop_codon'].reverse()
                for exon in rnaDict['exons']:
                    exon.reverse()

                rnaDict['exons'] = sorted(rnaDict['exons'], key=lambda exon : exon[0], reverse=True) 

            else:
                raise Exception("invalid strand : " + rnaDict['strand'])

    return annotDict

def getHgvsInfo(pos, var_g_p, txInfoList, chrAc, hdpConnexion):
    txExonsList = hdpConnexion.get_tx_exons(txInfoList[3], chrAc, 'splign')
    variantmapper = hgvs.variantmapper.EasyVariantMapper(hdpConnexion, primary_assembly='GRCh37')

    minDist = None
    for exonList in txExonsList:
        dist = min(abs(pos - exonList[8]), abs(pos - exonList[9]))
        
        if minDist == None: # initialisation
            minDist = dist
            exonNum = int(exonList[5]) + 1
        else:
            if dist < minDist:
                minDist = dist
                exonNum = int(exonList[5]) + 1

    var_c_p = variantmapper.g_to_c(var_g_p, txInfoList[3])

    return var_c_p, exonNum, minDist

def parsePileupLine(pileupLine):
    chrom, pos, ref, depth, readBase, qual = pileupLine.rstrip('\n').split('\t')
    
    ref = ref.upper()

    depth = 0
    
    countA = 0
    countT = 0
    countG = 0
    countC = 0
    countN = 0
    counta = 0
    countt = 0
    countg = 0
    countc = 0
    countn = 0
    insCountDict = dict()
    insString = '0'
    delCountDict = dict()
    delString = '0'

    # clean read base
    readBase = re.sub('\^.', '', readBase)
    readBase = re.sub('\$', '', readBase)

    # find insertion and deletion
    insMatch = re.findall('(\+[ATCGNatcgn0-9]+)', readBase)
    delMatch = re.findall('(\-[ATCGNatcgn0-9]+)', readBase)

    # handle insertion
    for ins in insMatch:
        depth += 1
        insLen = int(re.findall('\+([0-9]+)', ins)[0])
        insSeq = re.findall('\+' + str(insLen) + '([ATCGNatcgn]{' + str(insLen) + '})', ins)[0]
        # clean readbase
        readBase = re.sub('\+' + str(insLen) + insSeq, '', readBase)

        # add count to count dict
        if insSeq in insCountDict.keys():
            insCountDict[insSeq][0] += 1
        elif insSeq.upper() in insCountDict.keys():
            insCountDict[insSeq.upper()][1] += 1
        elif insSeq.isupper():
            insCountDict[insSeq] = [1, 0]
        else:
            insCountDict[insSeq.upper()] = [0, 1]


    # handle deletion
    for deletion in delMatch:
        depth += 1
        delLen = int(re.findall('\-([0-9]+)', deletion)[0])
        delSeq = re.findall('\-' + str(delLen) + '([ATCGNatcgn]{' + str(delLen) + '})', deletion)[0]
        # clean readbase
        readBase = re.sub('\-' + str(delLen) + delSeq, '', readBase)

        # add count to count dict
        if delSeq in delCountDict.keys():
            delCountDict[delSeq][0] += 1
        elif delSeq.upper() in delCountDict.keys():
            delCountDict[delSeq.upper()][1] += 1
        elif delSeq.isupper():
            delCountDict[delSeq] = [1, 0]
        else:
            delCountDict[delSeq.upper()] = [0, 1]

    # count snv
    for base in readBase:
        depth += 1
        # get proper base for reference supporting reads
        if base == '.':
            base = ref
        elif base == ',':
            base = ref.lower()

        if base == 'A': 
            countA += 1
        elif base =='a':
            counta += 1
        elif base == 'T':
            countT += 1
        elif base == 't':
            countt += 1
        elif base == 'C':
            countC += 1
        elif base == 'c':
            countc += 1
        elif base == 'G':
            countG += 1
        elif base == 'g':
            countg += 1
        elif base == 'N' or base == 'n':
            countN += 1
        elif base == '*':
            if base in delCountDict.keys():
                delCountDict[base][0] += 1
            else:
                delCountDict[base] = [1]
        else:
            raise Exception("invalid base in read base field : " + base + '\nline :' + line)

    pileupDict = {
        'chrom' : chrom, 'pos': pos,
        'depth': depth, 'ref' : ref,
        'A' : [countA, counta], 
        'T' : [countT, countt], 
        'C' : [countC, countc], 
        'G' : [countG, countg],
        'countN' : countN, 
        'ins' : insCountDict, 'del' : delCountDict
    }
    return pileupDict

def getCount(chrom, pos, refSeq, altSeq, mapQThreshold, baseQThreshold, refFastaFileName, bamFileName):

    """Get allelic ratio, depth and bases counts

    :param chrom: the chromosome where the variat is
    :param pos: the position of the variant
    :param altSeq: the sequence of the variant
    :param mapQThreshold:
    :param baseQThreshold:
    :param refFastaFileName:
    :param bamFileName:
    :return freq, depth, varCount, varCountF, varCountR: 
    """

    try:
        lenRef = len(refSeq)
        lenAlt = len(altSeq)

        tempPileupFile = tempfile.NamedTemporaryFile(mode="a+b", delete = False)
        tempFileName = tempPileupFile.name
        cmd = ' '.join(["samtools", "mpileup", "-AB", "-f", refFastaFileName, 
            "-q", str(mapQThreshold), "-Q", str(baseQThreshold), "-d 1000000",
            "-r", chrom + ':' + str(pos) + '-' + str(pos), bamFileName]
        )
        sp = subprocess.Popen(cmd, stdout = tempPileupFile, stderr=STDERR_FILE, shell=True)
        sp.wait()
        tempPileupFile.close()
        pileupFile = open(tempFileName, 'r')

        pileupLine = pileupFile.readline()
        pLineDict = parsePileupLine(pileupLine)

        # initialize counters
        varCount = 0
        varCountF = 0
        varCountR = 0

        # handle snv
        if lenRef == lenAlt:
            varCount = sum(pLineDict[altSeq])
            varCountF = pLineDict[altSeq][0]
            varCountR = pLineDict[altSeq][1]

        # handle del
        elif lenRef > lenAlt : 
            diff =  lenRef - lenAlt
            varCount = 0
            for delSeq, countList in zip(pLineDict['del'].keys(), pLineDict['del'].values()):
                if len(delSeq) >= diff and delSeq != '*':
                    varCount += sum(countList)
                    varCountF += countList[0]
                    varCountR += countList[1]

        # handle ins
        elif lenRef < lenAlt: 
            diff = lenAlt - lenRef
            varCount = 0
            for insSeq, countList in zip(pLineDict['ins'].keys(), pLineDict['ins'].values()):
                if insSeq[:diff].upper() == altSeq[lenRef:].upper():
                    varCount += sum(countList)
                    varCountF += countList[0]
                    varCountR += countList[1]

        if pLineDict['depth'] != 0:
            freq = '%.2f' % ((float(varCount)/float(pLineDict['depth'])) * 100)
        else:
            freq = "NA"

    finally:
        if tempFileName:
            pileupFile.close()
            os.remove(tempFileName)

    # get strand biases
    try:
        refStrandBias = float(max(pLineDict[refSeq[0]]))/float(sum(pLineDict[refSeq[0]]))
    except ZeroDivisionError:
        refStrandBias = 0

    try:
        varStrandBias = float(max(varCountF, varCountR))/float(varCountF + varCountR)
    except ZeroDivisionError:
        varStrandBias = 0
    
    return freq, pLineDict['depth'], sum(pLineDict[refSeq[0]]), pLineDict[refSeq[0]][0], pLineDict[refSeq[0]][1], '%.2f' % refStrandBias, varCount, varCountF, varCountR, '%.2f' % varStrandBias, pLineDict

def getCountEachBase(pLineDict):

    outList = list()
    if pLineDict != None:
        for base in 'ATCG':
            countFw = pLineDict[base][0]
            countRv = pLineDict[base][1]
            try:
                strandBias = float(max(countFw, countRv))/float(countFw + countRv)
            except ZeroDivisionError:
                strandBias = 0

            if pLineDict['depth'] == 0:
                allelicRatio = 0
            else:
                allelicRatio = (float(countFw + countRv)/float(pLineDict['depth'])) * 100

            outList.extend([countFw, countRv, '%.2f' % allelicRatio, '%.2f' % strandBias])
    else:
        logging.warn("No reads at position : counts and allelic ratio will be set to 0")
        for base in 'ATCG':
            outList.extend([0, 0, '%.2f' % 0, '%.2f' % 0])

    return outList

def getStatus(freq, threshold):
    if freq == "NA":
        status = "NA"
    elif float(freq) >= threshold:
        status = "Homozygous"
    else:
        status = "Heterozygous"

    return status

def getGranthamScore(refAA, altAA, prefered_nm, granthamDict):
    ## get and add grantham score
    if prefered_nm != "NA" and altAA and refAA in granthamDict.keys() and altAA in granthamDict[refAA]:
        return granthamDict[refAA][altAA]
    else:
        return "NA"

def getSpliceDist(chrom, pos, prefered_nm, annotDict):

        if prefered_nm == "NA":
            return "NA"
        else:
            spliceDistList = list()
            exonList = annotDict[chrom][prefered_nm]['exons']
            for exon in exonList:
                spliceDistList.append(abs(pos - int(exon[0])))
                spliceDistList.append(abs(pos - int(exon[1])))

            if len(spliceDistList) != 0:
                return min(spliceDistList)
            else:
                return "NA"

def median(numList):
    numList = sorted(numList)
    if len(numList) < 1:
        return None
    if len(numList) %2 == 1:
        return float(numList[((len(numList)+1)/2) - 1])
    else:
        return float(sum(numList[(len(numList)/2) - 1:(len(numList)/2) + 1]))/2.0

def getDistFromReadsStartEnd(chrom, pos, BamFile):
    pileupIt = BamFile.pileup(chrom, pos, pos + 1)
    endDistList = list()
    startDistList = list()

    for PileupColumn in pileupIt:
        for PileupRead in PileupColumn.pileups:
            if not PileupRead.alignment.is_secondary:
                endDistList.append(abs(pos - PileupRead.alignment.aend))
                startDistList.append(abs(pos - PileupRead.alignment.pos))

    return max(startDistList), median(startDistList), min(startDistList), max(endDistList), median(endDistList), min(endDistList)

#def getVarPositionInReads(chrom, pos, BamFile):
#    pileupIt = BamFile.pileup(chrom, pos, pos + 1)
#    locationList = list()
#
#    for PileupColumn in pileupIt:
#        for PileupRead in PileupColumn.pileups:
#            if not PileupRead.alignment.is_secondary and not PileupRead.alignment.is_unmapped and not PileupRead.alignment.mate_is_unmapped and PileupRead.alignment.mate_is_reverse:

def getMesScores(chrom, pos, ref, alt, exonNumber, annotDict, refFasta, prefered_nm, mesScriptDir):

    refMesScore = 'NA'
    altMesScore = 'NA'
    deltaMesScore = 'NA'
    chromAnnotDict = annotDict[chrom]
    if prefered_nm != 'NA':
        try: 
            nmDict = chromAnnotDict[prefered_nm]
        except KeyError:
            return "NA", "NA", "NA"
        # get minimum distance from splicing acceptor (index = 0)/donor (index = 1) site
        refSeq = None
        altSeq = None
        lenExons = len(nmDict['exons'])
        minDist = None
        index = None
        for i in xrange(0, lenExons):
            exon = nmDict['exons'][i]
            threePDist = abs(pos - exon[0])
            fivePDist = abs(pos - exon[1])

            if minDist == None:
                minDist = min(threePDist, fivePDist)
                exonNumber = i + 1
            else:
                if threePDist < minDist:
                    index = 0
                    minDist = threePDist
                    exonNumber = i + 1
                if fivePDist < minDist:
                    index = 1
                    minDist = fivePDist
                    exonNumber = i + 1


        exon = [elt + 1 for elt in nmDict['exons'][exonNumber - 1]]
        threePDist = abs(pos - exon[0])
        fivePDist = abs(pos - exon[1])

        # variant impacting acceptor sites
        if threePDist < fivePDist:
            if nmDict['strand'] == '+':
                dist = pos - exon[0] - 1
                refSeq = refFasta.fetch(chrom, exon[0] - 21, exon[0] + 2)

                if -21 <= dist < 2:
                    # check ref base
                    if ref[0].upper() != refSeq[21+dist].upper():
                        raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[21+dist]))

                    # get alternative sequence
                    if len(ref) == len(alt): # snv
                        altSeq = refSeq[:21 + dist] + alt + refSeq[21 + dist + 1:]
                    elif len(alt) > 1: # insertion
                        insSize = len(alt) - 1
                        altSeq = (refSeq[:21 + dist] + alt + refSeq[21 + dist + 1:])[-23:]
                    elif len(ref) > 1: # deletion
                        delSize = len(ref) - 1
                        if dist < 0:
                            altSeq = refFasta.fetch(chrom, exon[0] - 21 - delSize, exon[0] + 2)
                            altSeq = altSeq[delSize:21 + delSize + dist + 1] + alt + altSeq[21 + dist + delSize + 2:]
                        else:
                            altSeq = refFasta.fetch(chrom, exon[0] - 21, exon[0] + 4 + delSize)
                            altSeq = altSeq[:21 + dist] + alt + altSeq[21 + dist + 1 + delSize:]

                    cmd = "echo " + refSeq + " | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    refMesScore = float(process.stdout.read().rstrip())

                    cmd = "echo " + altSeq + " | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    altMesScore = float(process.stdout.read().rstrip())

            elif nmDict['strand'] == '-':
                dist = pos - exon[0] - 1
                refSeq = refFasta.fetch(chrom, exon[0] - 3, exon[0] + 20)
                #if -3 < dist <= 20:
                if 20 > dist >= -3:
                    if ref[0].upper() != refSeq[3 + dist].upper():
                        altSeq = refSeq[:3 + dist] + alt + refSeq[(3 + dist) + 1:]
                        raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[3+dist]))

                    # get alternative sequence
                    if len(ref) == len(alt): # snv
                        altSeq = refSeq[:3 + dist] + alt + refSeq[3 + dist + 1:]
                    elif len(alt) > 1: # insertion
                        insSize = len(alt) - 1
                        altSeq = (refSeq[:3 + dist] + alt + refSeq[3 + dist + 1:-insSize])[:23]
                    elif len(ref) > 1: # deletion
                        delSize = len(ref) - 1
                        if dist < 0:
                            altSeq = refFasta.fetch(chrom, exon[0] - 3 - delSize, exon[0] + 20)
                            altSeq = altSeq[delSize:3 + delSize + dist + 1] + alt + altSeq[3 + dist + delSize + 2:]
                        else:
                            altSeq = refFasta.fetch(chrom, exon[0] - 3, exon[0] + 20 + delSize)
                            altSeq = altSeq[:3 + dist] + alt + altSeq[3 + dist + 1 + delSize:]

                    # reverse and translate sequences
                    refSeq = refSeq.translate(TRANSLATION_TAB)[::-1]
                    altSeq = altSeq.translate(TRANSLATION_TAB)[::-1]

                    cmd = "echo " + refSeq + " | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    refMesScore = float(process.stdout.read().rstrip())

                    cmd = "echo " + altSeq + " | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    altMesScore = float(process.stdout.read().rstrip())

        # variant impacting donor sites
        else:
            if nmDict['strand'] == '-':
                dist = pos - exon[1] - 1
                refSeq = refFasta.fetch(chrom, exon[1] - 7, exon[1] + 2)
                if -7 <= dist < 2:
                    # check ref base
                    if ref[0].upper() != refSeq[7 + dist].upper():
                        raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[7+dist]))

                    # get alternative sequence
                    if len(ref) == len(alt): # snv
                        altSeq = refSeq[:7 + dist] + alt + refSeq[7 + dist + 1:]
                    elif len(alt) > 1: # insertion
                        insSize = len(alt) - 1
                        altSeq = (refSeq[:7 + dist] + alt + refSeq[7 + dist + 1:-insSize])[:23]
                    elif len(ref) > 1: # deletion
                        delSize = len(ref) - 1
                        if dist < 0:
                            altSeq = refFasta.fetch(chrom, exon[0] - 7 - delSize, exon[0] + 2)
                            altSeq = altSeq[delSize:7 + delSize + dist + 1] + alt + altSeq[7 + dist + delSize + 2:][-9:]
                        else:
                            altSeq = refFasta.fetch(chrom, exon[0] - 7, exon[0] + 2 + delSize)
                            altSeq = altSeq[:7 + dist] + alt + altSeq[7 + dist + 1 + delSize:]

                    refSeq = refSeq.translate(TRANSLATION_TAB)[::-1]
                    altSeq = altSeq.translate(TRANSLATION_TAB)[::-1]

                    cmd = "echo " + refSeq + " | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    refMesScore = float(process.stdout.read().rstrip())

                    cmd = "echo " + altSeq + " | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    altMesScore = float(process.stdout.read().rstrip())
                    
            elif nmDict['strand'] == '+':
                dist = pos - exon[1] - 1
                refSeq = refFasta.fetch(chrom, exon[1] - 3, exon[1] + 6)
                if -3 <= dist < 6:
                    # check ref base
                    if ref[0].upper() != refSeq[3 + dist].upper():
                        altSeq = refSeq[:3 + dist] + alt + refSeq[(3 + dist) + 1:]
                        raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[3 + dist]))

                    # get alternative sequence
                    if len(ref) == len(alt): # snv
                        altSeq = refSeq[:3 + dist] + alt + refSeq[3 + dist + 1:]
                    elif len(alt) > 1: # insertion
                        insSize = len(alt) - 1
                        altSeq = (refSeq[:3 + dist] + alt + refSeq[3 + dist + 1:-insSize])[:9]
                    elif len(ref) > 1: # deletion
                        delSize = len(ref) - 1
                        if dist < 0:
                            altSeq = refFasta.fetch(chrom, exon[1] - 3 - delSize, exon[1] + 6)
                            altSeq = altSeq[delSize:3 + delSize + dist + 1] + alt + altSeq[3 + dist + delSize + 2:]
                        else:
                            altSeq = refFasta.fetch(chrom, exon[1] - 3, exon[1] + 6 + delSize)
                            altSeq = altSeq[:3 + dist] + alt + altSeq[3 + dist + 1 + delSize:]

                    cmd = "echo " + refSeq + " | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    refMesScore = float(process.stdout.read().rstrip())

                    cmd = "echo " + altSeq + " | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2"
                    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
                    altMesScore = float(process.stdout.read().rstrip())

    if refMesScore < 0:
        refMesScore = 0
    if altMesScore < 0:
        altMesScore = 0

    if refMesScore != "NA" and altMesScore != "NA":
        if refMesScore == 0:
            logging.warning("reference maxentscan score is equal to 0 for position {0} on chormosome {1}. Delta mes cannot be processed, 'NA' will be reported".format(pos, chrom))
        else:
            deltaMesScore = altMesScore/refMesScore

    return refMesScore, altMesScore, deltaMesScore

def getMesSlidingWindow(chrom, ref, alt, pos, refFasta, mesScriptDir):

        # acceptor sites
        if len(ref) == 1 and len(alt) == 1:
            refSeq = refFasta.fetch(chrom, pos - 23, pos + 23)
            altSeq = ""
            swRange = xrange(0,23)
            if ref.upper() != refSeq[22].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[22]))
            else:
                altSeq = refSeq[:22] + alt + refSeq[23:]

        elif len(ref) > 1: # deletion
            delSize = len(ref) - len(alt)
            refSeq = refFasta.fetch(chrom, pos - 23, pos + 23)
            altSeq = refFasta.fetch(chrom, pos - 23, pos + 23 + delSize)
            swRange = xrange(0,23)
            if ref[0].upper() != refSeq[22].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[22]))
            else:
                altSeq = altSeq[:22] + alt + altSeq[23 + delSize:]

        elif len(alt) > 1: # insertion
            insSize = len(alt) - len(ref)
            refSeq = refFasta.fetch(chrom, pos - 23, pos + 23 + insSize)
            altSeq = refFasta.fetch(chrom, pos - 23, pos + 23)
            if ref.upper() != refSeq[22].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[22]))
            else:
                altSeq = altSeq[:22] + alt + altSeq[23:]

            swRange = xrange(0,23 + insSize)

        refSeqList = list()
        altSeqList = list()
        for i in swRange:
            refSeqList.append(refSeq[i:i+23])
            altSeqList.append(altSeq[i:i+23])
            
        cmd = "echo -e \"" + '\\n'.join(refSeqList) + "\" | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2 | tr '\\n' '\t'"
        process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
        refSWMescoreList = [float(x) for x in process.stdout.read().rstrip().split('\t')]
        refSWMescoreAcceptor = max(refSWMescoreList)

        cmd = "echo -e \"" + '\\n'.join(altSeqList) + "\" | perl " + mesScriptDir + "score3.pl /dev/stdin | cut -f2 | tr '\\n' '\t'"
        process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
        altSWMescoreList = [float(x) for x in process.stdout.read().rstrip().split('\t')]
        altSWMescoreAcceptor = max(altSWMescoreList)

        if refSWMescoreAcceptor < 0:
            refSWMescoreAcceptor = 0
        if altSWMescoreAcceptor < 0:
            altSWMescoreAcceptor = 0

        # donor sites
        if len(ref) == 1 and len(alt) == 1:
            refSeq = refFasta.fetch(chrom, pos - 9, pos + 9)
            altSeq = ""
            swRange = xrange(0,9)
            if ref.upper() != refSeq[8].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[8]))
            else:
                altSeq = refSeq[:8] + alt + refSeq[9:]

        elif len(ref) > 1: # deletion
            delSize = len(ref) - len(alt)
            refSeq = refFasta.fetch(chrom, pos - 9, pos + 9)
            altSeq = refFasta.fetch(chrom, pos - 9, pos + 9 + delSize)
            swRange = xrange(0,9)
            if ref[0].upper() != refSeq[8].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[8]))
            else:
                altSeq = altSeq[:8] + alt + altSeq[9 + delSize:]

        elif len(alt) > 1: # insertion
            insSize = len(alt) - len(ref)
            refSeq = refFasta.fetch(chrom, pos - 9, pos + 9 + insSize)
            altSeq = refFasta.fetch(chrom, pos - 9, pos + 9)
            if ref.upper() != refSeq[8].upper():
                raise Exception('input reference base does not match base in reference genome: position {0} input {1} ref {2}'.format(pos, ref, refSeq[8]))
            else:
                altSeq = altSeq[:8] + alt + altSeq[9:]

            swRange = xrange(0,9 + insSize)

        refSeqList = list()
        altSeqList = list()
        for i in swRange:
            refSeqList.append(refSeq[i:i+9])
            altSeqList.append(altSeq[i:i+9])
            
        cmd = "printf \"" + '\\n'.join(refSeqList) + "\" | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2 | tr '\\n' '\t'"
        process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
        refSWMescoreList = [float(x) for x in process.stdout.read().rstrip().split('\t')]
        refSWMescoreDonor = max(refSWMescoreList)

        cmd = "printf \"" + '\\n'.join(altSeqList) + "\" | perl " + mesScriptDir + "score5.pl /dev/stdin | cut -f2 | tr '\\n' '\t'"
        process = subprocess.Popen(cmd, stdout = subprocess.PIPE, cwd = mesScriptDir, shell=True)
        altSWMescoreList = [float(x) for x in process.stdout.read().rstrip().split('\t')]
        altSWMescoreDonor = max(altSWMescoreList)

        if refSWMescoreDonor < 0:
            refSWMescoreDonor = 0
        if altSWMescoreDonor < 0:
            altSWMescoreDonor = 0

        return refSWMescoreAcceptor, altSWMescoreAcceptor, refSWMescoreDonor, altSWMescoreDonor

def getESRScore(chrom, pos, refSeq, altSeq, esrScoreDict, refFasta):

        # get sequence
        refSubSeq = refFasta.fetch(chrom, pos - 6, pos +5)
        refLength = len(refSeq)
        altLength = len(altSeq)

        # handle indels
        indelSize = refLength - altLength
        if refLength != altLength:
            if indelSize < 0: # insertion
                # handle long insertions
                if indelSize <= -12:
                    altSubSeq = None
                else:
                    indelSubSeqList = list(refFasta.fetch(chrom, pos - 6, pos + 5 + indelSize))
                    altSubSeq = ''.join(indelSubSeqList[:5]) + altSeq + ''.join(indelSubSeqList[6:])
            else: # deletion
                indelSubSeqList = list(refFasta.fetch(chrom, pos - 6, pos + 5 + indelSize))
                altSubSeq = ''.join(indelSubSeqList[:6]) + ''.join(indelSubSeqList[6 + indelSize:])
        # handle snv
        else:
            altSubSeq = refSubSeq[:5] + altSeq + refSubSeq[6:]

        if altSubSeq != None:
            refScoreList = list()
            altScoreList = list()
            for i in xrange(0, len(refSubSeq) - 5):
                refHexamer = refSubSeq[i:i+6].upper()
                altHexamer = altSubSeq[i:i+6].upper()
                refESRscore = esrScoreDict[refHexamer]
                altESRscore = esrScoreDict[altHexamer]
                refScoreList.append(refESRscore)
                altScoreList.append(altESRscore)

            refScore = max(refScoreList)
            altScore = max(altScoreList)
            delta = altESRscore / refESRscore

        else:
            altScore = "NA"
            delta = "NA"

        return refScore, altScore, delta

