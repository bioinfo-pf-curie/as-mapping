#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   Script to map reads on parental genomes using Tophat

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
	exit
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

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_parental_tophat
mkdir -p ${BAM_OUT}
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_parental

# Checking input parameters for mapping
#  Indexes ?
if [[ -z ${B2_INDEX_GENO1} ]]; then B2_INDEX_GENO1=${INDEXES}/${ID_GENO1};fi
if [[ ! -e ${B2_INDEX_GENO1}.rev.2.bt2 ]]
then
	echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_GENO1}. Exit." 1>&2
    exit 1
fi
if [[ -z ${B2_INDEX_GENO2} ]]; then B2_INDEX_GENO2=${INDEXES}/${ID_GENO2};fi
if [[ ! -e ${B2_INDEX_GENO2}.rev.2.bt2 ]]
then
	echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_GENO2}. Exit." 1>&2
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

## Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    # - Mapping
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out_${ID_GENO1} ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_GENO1} ${FQ_READS_F}
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out_${ID_GENO2} ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_GENO2} ${FQ_READS_F}
    mv ${BAM_OUT}/tophat_out_${ID_GENO1}/accepted_hits.bam ${BAM_OUT}/tophat_out_${ID_GENO1}/${ID_GENO1}.bam
    mv ${BAM_OUT}/tophat_out_${ID_GENO2}/accepted_hits.bam ${BAM_OUT}/tophat_out_${ID_GENO2}/${ID_GENO2}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -p ${BAM_OUT}/tophat_out_${ID_GENO1}/${ID_GENO1}.bam -m ${BAM_OUT}/tophat_out_${ID_GENO2}/${ID_GENO2}.bam -o ${BAM_OUT} -n ${ID_OUTBAM}
else
    # Paired-end mapping
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out_${ID_GENO1} ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_GENO1} ${FQ_READS_F} ${FQ_READS_R}
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out_${ID_GENO2} ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_GENO2} ${FQ_READS_F} ${FQ_READS_R}
    mv ${BAM_OUT}/tophat_out_${ID_GENO1}/accepted_hits.bam ${BAM_OUT}/tophat_out_${ID_GENO1}/${ID_GENO1}.bam
    mv ${BAM_OUT}/tophat_out_${ID_GENO2}/accepted_hits.bam ${BAM_OUT}/tophat_out_${ID_GENO2}/${ID_GENO2}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -p ${BAM_OUT}/tophat_out_${ID_GENO1}/${ID_GENO1}.bam -m ${BAM_OUT}/tophat_out_${ID_GENO2}/${ID_GENO2}.bam -s 2 -o ${BAM_OUT} -n ${ID_OUTBAM}
fi


exit 0 
