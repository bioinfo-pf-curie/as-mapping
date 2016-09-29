#!/bin/bash
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
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_NMASK} ${FQ_READS_F}
else
    # Paired-end mapping
    ${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_NMASK} ${FQ_READS_F} ${FQ_READS_R}
fi


# Filter [TO BE DONE] using markAllelicStatus

exit 0
