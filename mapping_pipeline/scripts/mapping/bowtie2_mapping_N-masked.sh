#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   Script to perform alignment to a N-masked genome using Bowtie2
#   and filter for allele specific analysis

#### Parameters #### --------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-c) config=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: ERROR - unrecognized option $1" 1>&2; exit 1;; 
    (*)  break;;
    esac
    shift
done

if [[ -z $config ]]
then
    echo "$0: ERROR - you need to specify a config file. Exit." 1>&2
    exit 1
fi
source ${config}

#### Function #### ----------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------

## Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_N-masked_bowtie2
mkdir -p ${BAM_OUT}
ID_NMASK=N-masked_${ID_GENO1}_${ID_GENO2}
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_N-masked

## Checking input parameters for mapping
#- Indexes ?
if [[ -z ${B2_INDEX_NMASK} ]]; then B2_INDEX_NMASK=${INDEXES}/${ID_NMASK};fi
if [[ ! -e ${B2_INDEX_NMASK}.rev.2.bt2 ]]
then
    echo "$0: ERROR - Missing Bowtie2 indexes. Exit." 1>&2
    exit 1
fi
#- Reads ?
if [[ ! -e ${FQ_READS_F} ]]
then
    echo "$0: ERROR - Missing FASTQ reads. Exit." 1>&2
    exit 1
fi
SINGLE_END=true
if [[ -e ${FQ_READS_R} ]]; then SINGLE_END=false;fi
#- SNP file ?
if [[ -z ${SNP_FILE} ]]; then SNP_FILE=${OUT_DIR}/all_SNPs_${PATERNAL}_${MATERNAL}.txt.gz;fi
if [[ ! -e ${SNP_FILE} ]]
then
    echo "$0: WARNING - Missing SNP file for SNPsplit." 1>&2
fi

## Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_NMASK} -U ${FQ_READS_F} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    ${SNPSPLIT} --samtools_path ${SAMTOOLS_DIR} --snp_file ${SNP_FILE} ${BAM_OUT}/${ID_OUTBAM}.bam 
else
    # Paired-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_NMASK} -1 ${FQ_READS_F} -2 ${FQ_READS_R} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    ${SNPSPLIT} --paired --samtools_path ${SAMTOOLS_DIR} --snp_file ${SNP_FILE} ${BAM_OUT}/${ID_OUTBAM}.bam 
fi

# Removing files
rm ${BAM_OUT}/${ID_OUTBAM}.SNPsplit_sort.txt ${BAM_OUT}/${ID_OUTBAM}.genome*.bam ${BAM_OUT}/${ID_OUTBAM}.unassigned.bam ${BAM_OUT}/${ID_OUTBAM}.sortedByName.bam

exit 0
