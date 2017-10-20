#! /usr/bin/env python
# -*- coding: utf8 -*-

# This file is part of ClinTools
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

import re
import os
import pysam
import optparse
import logging

import clinTools

import hgvs
import hgvs.variantmapper
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.validator

# create usefuls objects for hgvs module use
# hgvs position parser
hgvsparser = hgvs.parser.Parser()
# create a variant validator
#vr = hgvs.validator.Validator(hdp=hdp)
# tips : to search availables nm in database for a gene use:
# hdp.get_tx_for_gene('nom_de_gene')

# lists of fields
fieldList = [
    "Barcode", 
    # rename add 'Gene.refGene'
    "Gene",
    # first 4 fields from annovar file [1:5]
    "Chr", "Start", "End", 
    # ref and alt get form vcf part of annovar file
    "Ref", "Alt", 
    # Fields corresponding to 'Func.refGene' and 'ExonicFunc.refGene' in annovar file
    "Location", "Type", 
    # other fields
    "NM", "cDNAchange", 
    # "AAchange", 
    "AAchange",
    "all_hgvs", "Status", "Depth", "Count_ref", "Count_ref_F", "Count_ref_R", "Ref_strand_bias",
    "Count_alt", "Count_alt_F", "Count_alt_R", "Alt_strand_bias", "Allelic_ratio_%", 
    "Count_A_F", "Count_A_R", "A_allelic_ratio", "A_strand_bias",
    "Count_T_F", "Count_T_R", "T_allelic_ratio", "T_strand_bias",
    "Count_C_F", "Count_C_R", "C_allelic_ratio", "C_strand_bias",
    "Count_G_F", "Count_G_R", "G_allelic_ratio", "G_strand_bias",
    # grantham
    "Grantham", 
    "can_splice_dist", "MES_ref", "MES_alt", "MES_delta",
    "exon", 
    "MES_cryptic_acceptor_wt", "MES_cryptic_acceptor_mut","MES_cryptic_donor_wt", "MES_cryptic_donor_mut",
    "Bases_around", 
    "max_distance_from_reads_start", "median_distance_from_reads_start", "min_distance_from_reads_start", 
    "max_distance_from_reads_end", "median_distance_from_reads_end", "min_distance_from_reads_end"
]

annovarMandatoryFields = ["Start", "End", "Chr", "Otherinfo", "Ref", "Alt", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene"]


## Get options and read useful conf files
parser = optparse.OptionParser(description="Generate an annotation report for a set of table from annovar",
    usage="usage: %prog [options] input_conf_file prefered_nm_file_conf annotation_gtf reference_fasta_file chr_accesion_file"
)
parser.add_option("-d", "--hgvs-max-dist", type=int, default=200,
    help="Maximum distance from start of first exon or end of last exon of an nm to report an hgvs position (default : 200)"
)
parser.add_option("-b", "--uta-database", type=str, default=None,
    help="The url of uta database to use (ex : postgresql://anonymous:anonymous@uta.biocommons.org/uta_dev/uta_20150704). If not set the default database is used."
)
parser.add_option("-q", "--mapQ", default=0,
    help="-q option for mpileup for frequency processing and base counting. mapQ threshold (default : 0)"
)
parser.add_option("-Q", "--BAQ", default=0,
    help="-Q option for mpileup for frequency processing and base counting. base quality threshold (default : 0)"
)
(options, args) = parser.parse_args()

if len(args) < 5:
    parser.print_help()
    exit(1)

clinTools.checkSamtoolsVersion('1.1')

inputConfFileName = args[0]
preferedNmFile = args[1]
gtfFileName = args[2]
refFastaFileName = args[3]
chrAccesionFile = args[4]

mesScriptDir = os.path.dirname(os.path.realpath(__file__)) + "/maxentscan/"
scoreTablePath = os.path.dirname(os.path.realpath(__file__)) + "/score_tables/"

granthamFileName = scoreTablePath + "grantham.tsv"

hgvsMaxDist = options.hgvs_max_dist
REF_FASTA = pysam.Fastafile(refFastaFileName)

# connect to database
if options.uta_database:
    hdp = hgvs.dataproviders.uta.connect(db_url=options.uta_database)
else:
    hdp = hgvs.dataproviders.uta.connect()

# create variantmapper to convert between g. c. p.
variantmapper = hgvs.variantmapper.EasyVariantMapper(hdp, primary_assembly='GRCh37')

# create annotation dictionnary from gtf file (use for maxentscan score and hgvs)
annotDict = clinTools.parseChrGtf(gtfFileName)

# create prefered nm list
preferedNmList = list()
with open(preferedNmFile, 'r') as pNmF:
    lines = pNmF.readlines()
    for line in lines:
        preferedNmList.append(line.rstrip().split('\t')[0])
        
# get grantham dict from file
granthamDict = dict()
with open(granthamFileName, 'r') as granthamFile:
    gHeader = granthamFile.readline().rstrip().split('\t')
    gLines = granthamFile.readlines()
    for line in gLines:
        lineList = line.rstrip().split()
        granthamDict[lineList[0]] = dict()
        for i in xrange(1, len(lineList)):
            granthamDict[lineList[0]][gHeader[i]] = lineList[i]

# get chromosome accession dict from file
# this dict is useful to format hgvs position for hgvs module
chrAcDict = dict()
with open(chrAccesionFile, 'r') as chrAcFile:
    for line in chrAcFile:
        lineList = line.rstrip().split('\t')
        chrAcDict[lineList[0]] = lineList[1]

# generate header line and get annovar field from first file
annovarFieldList = list()
with open(inputConfFileName, 'r') as inputConfFile:
    headerList = fieldList[:9]
    annovarFile = open(inputConfFile.readlines()[0].rstrip().split('\t')[1] , 'r')
    lineList = annovarFile.readlines()[0].rstrip().split('\t')
    for field in lineList:
        if not field in annovarMandatoryFields:
            headerList.append(field)
            annovarFieldList.append(field)
    headerList.extend(fieldList[9:])

# print header
print '\t'.join(headerList)

## Begin to do something on annovar files
inputConfFile = open(inputConfFileName, 'r')
for confLine in inputConfFile:
    barcode, annovarFile, bamFile = confLine.rstrip().split('\t')

    annotatedFile = open(annovarFile, 'r')
    BamFile = pysam.AlignmentFile(bamFile, 'rb')
    lines = annotatedFile.readlines()
    # get field order in annovar file
    aFieldList = lines[0].rstrip().split('\t')

    # check if all mandatory annovar fields are present
    for field in annovarMandatoryFields:
        if not field in aFieldList:
            raise Exception('Mandatory field "{0}" is not in file "{1}" barcode "{2}"'.format(field, annovarFile, barcode))

    for line in lines[1:]:
        lineList = line.rstrip().split('\t')

        pos = int(lineList[aFieldList.index("Otherinfo") + 1])
        end = int(lineList[aFieldList.index("End")])
        chrom = lineList[aFieldList.index("Chr")]
        ref = lineList[aFieldList.index("Otherinfo") + 3]
        alt = lineList[aFieldList.index("Otherinfo") + 4]
        refLength = len(ref)
        altLength = len(alt)

        availableTxList = hdp.get_tx_for_region(chrAcDict[chrom], 'splign', pos, pos)
        prefered_nm = "NA"
        prefered_nm_v = "NA"
        lenAlt = len(alt)
        lenRef = len(ref)

        # select nm
        for txList in availableTxList:
            for nm in preferedNmList:
                if re.search(nm, txList[0]):
                    prefered_nm_v = txList[0]
                    prefered_nm = nm

        if prefered_nm == "NA":
            if lineList[aFieldList.index("AAChange.refGene")] != 'NA':
                prefered_nm = lineList[aFieldList.index("AAChange.refGene")].split(',')[0].split(':')[1]
                for txList in availableTxList:
                    if re.search(prefered_nm, txList[0]):
                        prefered_nm_v = txList[0]

        ## get hgvs from hgvs module for proper p. conversion
        # generate g. for hgvs module
        if lenAlt > lenRef: # insertion
            g_point = '{ac}:g.{pos}ins{alt}'.format(ac=chrAcDict[chrom], pos=pos + 1, alt=alt[1:])
        elif lenRef > lenAlt: # deletion
            g_point = '{ac}:g.{pos}del{ref}'.format(ac=chrAcDict[chrom], pos=pos + 1, ref=ref[1:])
        else:
            g_point = '{ac}:g.{pos}{ref}>{alt}'.format(ac=chrAcDict[chrom], pos=pos, ref=ref, alt=alt)

        var_g_p = hgvsparser.parse_hgvs_variant(g_point)

        # get hgvs for prefered nm
        spliceDist = "NA"
        if prefered_nm_v == "NA":
            c_point = "NA"
            p_point = "NA"
            exonNumber = "NA"

        else:
            try:
                txInfoList = hdp.get_tx_info(prefered_nm_v, chrAcDict[chrom],'splign')
                var_c_p, exonNumber, spliceDist = clinTools.getHgvsInfo(pos, var_g_p, txInfoList, chrAcDict[chrom], hdp)
                var_p_p = variantmapper.c_to_p(var_c_p)
                c_point = str(var_c_p).split(':')[-1]
                p_point = str(var_p_p).split(':')[-1]
            except hgvs.exceptions.HGVSInvalidIntervalError:
                c_point = "NA"
                p_point = "NA"
                exonNumber = "NA"

        # report c. and p. for all nm available at position
        allHgvsList = list()
        for txList in availableTxList:
            try:
                txInfoList = hdp.get_tx_info(txList[0], chrAcDict[chrom],'splign')
                var_c_p, exonNumber, spliceDist = clinTools.getHgvsInfo(pos, var_g_p, txInfoList, chrAcDict[chrom], hdp)
                var_p_p = variantmapper.c_to_p(var_c_p)
                allHgvsList.append(str(var_c_p) + "," + str(var_p_p))
            except hgvs.exceptions.HGVSInvalidIntervalError:
                pass

        if len(allHgvsList) == 0:
            allHgvsList.append("NA")
            
        if p_point != "NA" and var_p_p.posedit.edit != '=' and var_p_p.posedit.edit != '?':
            refAA = var_p_p.posedit.pos.start.aa
            altAA = var_p_p.posedit.edit.alt
        else:
            refAA = "NA"
            altAA = "NA"

        # calculate coverage and allelic ratio from bam 
        try:
            freq, depth, refCount, refCountF, refCountR, refStrandBias, varCount, varCountF, varCountR, varStrandBias, pLineDict  = clinTools.getCount(
                chrom, pos, ref, alt, options.mapQ, options.BAQ, refFastaFileName, bamFile
            )
        except ValueError as e:
            freq, depth, refCount, refCountF, refCountR, refStrandBias, varCount, varCountF, varCountR, varStrandBias = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            pLineDict = None
            logging.warn("variant not found in bam. 0 count and allelic ratio will be reported")

        # status Heterozygote or Homozygotes
        status = clinTools.getStatus(freq, 75)
        # add grantham scores
        granthamScore = clinTools.getGranthamScore(refAA, altAA, prefered_nm, granthamDict)

        ## maxentscan for canonic splicing site
        refMesScore, altMesScore, deltaMesScore = clinTools.getMesScores(chrom, pos, ref, alt, exonNumber, annotDict, REF_FASTA, prefered_nm, mesScriptDir)
        ## maxentscan sliding windows --> report max scores
        try:
            refSWMescoreAcceptor, altSWMescoreAcceptor, refSWMescoreDonor, altSWMescoreDonor = clinTools.getMesSlidingWindow(chrom, ref, alt, pos, REF_FASTA, mesScriptDir)
        except ValueError as e:
            refSWMescoreAcceptor, altSWMescoreAcceptor, refSWMescoreDonor, altSWMescoreDonor = 0, 0, 0, 0
            logging.warn("no mes on sliding windows available")

        # get base around
        baseAround = REF_FASTA.fetch(chrom, pos - 15, pos + 15)
        # calculate min and max distances between variant position and reads start / ends
        try:
            maxDistFromStart, medianDistFromStart, minDistFromStart, maxDistFromEnd, medianDistFromEnd, minDistFromEnd = clinTools.getDistFromReadsStartEnd(chrom, pos, BamFile)
        except ValueError:
            maxDistFromStart, medianDistFromStart, minDistFromStart, maxDistFromEnd, medianDistFromEnd, minDistFromEnd = 0, 0, 0, 0, 0, 0
            logging.warn("No distance from read can be calculated 0 will be reported for each values")


        ## add  "Gene" "Chr", "Start", "End", "Ref", "Alt", "Location", "Type"
        outLineList = [barcode, lineList[aFieldList.index('Gene.refGene')], chrom, pos, end, ref, alt, lineList[aFieldList.index('Func.refGene')], lineList[aFieldList.index('ExonicFunc.refGene')]]

        # add field from annovar 
        for field in annovarFieldList:
            try:
                value = lineList[aFieldList.index(field)]
            except ValueError:
                raise Exception('Field "{0}" is not in file "{1}" barcode "{2}". All annovar files must contain the same fields'.format(field, annovarFile, barcode))
            outLineList.append(value)

        outLineList.extend([prefered_nm_v, c_point, p_point])
        ## add all_hgvs
        outLineList.append(';'.join(allHgvsList))
        # add status, depth of coverages, freq and strand bias for called variant
        outLineList.extend([status, depth, refCount, refCountF, refCountR, refStrandBias, varCount, varCountF, varCountR, varStrandBias, freq])
        # add count, strand bias and allelic ratio for each base
        outLineList.extend(clinTools.getCountEachBase(pLineDict))
        outLineList.append(granthamScore)
        ## add canonic splice distance
        outLineList.append(spliceDist)

        outLineList.extend([refMesScore, altMesScore, deltaMesScore])

        # add exon number 
        outLineList.append(exonNumber)
        outLineList.extend([refSWMescoreAcceptor, altSWMescoreAcceptor, refSWMescoreDonor, altSWMescoreDonor])

        outLineList.append(baseAround)
        outLineList.extend([maxDistFromStart, medianDistFromStart, minDistFromStart, maxDistFromEnd, medianDistFromEnd, minDistFromEnd])

        # print result line
        print '\t'.join([str(elt) for elt in outLineList])

inputConfFile.close()
annotatedFile.close()
