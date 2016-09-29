#!/bin/bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to perform alignment to a diploi genome

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

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_diploid_bowtie2
mkdir -p ${BAM_OUT}
ID_DIP=${ID_GENO1}_${ID_GENO2}
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_diploid

## Checking input parameters for mapping
#- Indexes ?
if [[ -z ${B2_INDEX_DIP} ]]; then B2_INDEX_DIP=${INDEXES}/${ID_DIP};fi
if [[ ! -e ${B2_INDEX_DIP}.rev.2.bt2l ]]
then
    echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_DIP}. Exit." 1>&2
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

## Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    # - Mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_DIP} -U ${FQ_READS_F} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -d ${BAM_OUT}/${ID_OUTBAM}.bam -o ${BAM_OUT} -n ${ID_OUTBAM}
else
    # Paired-end mapping
    # - Mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_DIP} -1 ${FQ_READS_F} -2 ${FQ_READS_R} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -d ${BAM_OUT}/${ID_OUTBAM}.bam -s 2 -o ${BAM_OUT} -n ${ID_OUTBAM}
fi


exit 0
