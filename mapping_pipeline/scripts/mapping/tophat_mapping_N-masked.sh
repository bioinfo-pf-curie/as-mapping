#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   Script to perform alignment to a N-masked genome using Tophat
#   and filter for allele specific analysis

#### Parameters #### --------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-c) config=$2; shift;;
    (-f) FQ_READS_F=$2; shift;;
    (-r) FQ_READS_R=$2; shift;;
    (-o) OUT_DIR=$2; shift;;
    (-n) OUT_NAME=$2; shift;;
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
    echo -e "-f"" <Forward reads (paired-end) or reads (single-end)>"
    echo -e "-r"" <[Paired-end ONLY] reverse reads>"
    echo -e "-o"" <Output directory for the alignment files>"
    echo -e "-n"" <Output name for the alignment files>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------

## Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_N-masked_tophat
mkdir -p ${BAM_OUT}
ID_NMASK=N-masked_${ID_GENO1}_${ID_GENO2}
# Name for the output BAM file
ID_OUTBAM=${OUT_NAME}_N-masked

## Checking input parameters for mapping
#  Indexes ?
if [[ -z ${B2_INDEX_NMASK} ]]; then B2_INDEX_NMASK=${INDEXES}/${ID_NMASK}; fi
if [[ ! -e ${B2_INDEX_NMASK}.rev.2.bt2 ]]
then
    echo "$0: ERROR - Missing Bowtie2 indexes. Exit." 1>&2
    echo "$0: ERROR - HINT: Generate genomes and indexes using build_genomes_indexes.sh." 1>&2
    exit 1
fi
#  Reads ?
if [[ ! -e ${FQ_READS_F} ]]
then
    echo "$0: ERROR - Missing FASTQ reads. Exit." 1>&2
    exit 1
fi
SINGLE_END=true
if [[ -e ${FQ_READS_R} ]]; then SINGLE_END=false;fi

#  Tophat annotation ?
if [[ -e ${TOPHAT_TRANSC_INDEX} ]]
then
    TOPHAT_OPT="${TOPHAT_OPT} --transcriptome-index ${TOPHAT_TRANSC_INDEX}"
else
    if [[ -e ${TOPHAT_GTF} ]]
    then
        TOPHAT_OPT="${TOPHAT_OPT} --GTF ${TOPHAT_GTF}"
    else
        echo "$0: ERROR - Missing transcript annotation. Exit." 1>&2
        exit 1
    fi  
fi
#- SNP file ?
if [[ -z ${SNP_FILE} ]]; then SNP_FILE=${FASTA_OUT}/all_SNPs_${PATERNAL}_${MATERNAL}.txt.gz;fi
if [[ ! -e ${SNP_FILE} ]]
then
    echo "$0: WARNING - Missing SNP file for SNPsplit." 1>&2
fi

## Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_NMASK} ${FQ_READS_F}
    mv ${BAM_OUT}/tophat_out/accepted_hits.bam ${BAM_OUT}/${ID_OUTBAM}.bam
    ${SNPSPLIT} --samtools_path ${SAMTOOLS_DIR} --snp_file ${SNP_FILE} ${BAM_OUT}/${ID_OUTBAM}.bam
else
    # Paired-end mapping
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_NMASK} ${FQ_READS_F} ${FQ_READS_R}
    mv ${BAM_OUT}/tophat_out/accepted_hits.bam ${BAM_OUT}/${ID_OUTBAM}.bam
    ${SNPSPLIT} --paired --samtools_path ${SAMTOOLS_DIR} --snp_file ${SNP_FILE} ${BAM_OUT}/${ID_OUTBAM}.bam
fi

# Removing files
rm ${BAM_OUT}/${ID_OUTBAM}.SNPsplit_sort.txt ${BAM_OUT}/${ID_OUTBAM}.genome*.bam ${BAM_OUT}/${ID_OUTBAM}.unassigned.bam ${BAM_OUT}/${ID_OUTBAM}.sortedByName.bam

exit 0
