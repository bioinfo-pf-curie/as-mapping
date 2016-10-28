#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to perform alignment to parental genomes

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
source ${PIPELINE_PATH}/includes/path_fct.inc


#### Main #### --------------------------------------------------------------------------

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_parental_bowtie2
mkdir -p ${BAM_OUT}
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_parental

## Checking input parameters for mapping
#- One parent is reference ?
if [[ ${ID_GENO1} == 'REF' ]]
then
    ID_REF=$(basename ${REF_GENO%.fa*})
    if [[ -z ${B2_INDEX_REF} ]]
    then
        B2_INDEX_GENO1=${INDEXES}/${ID_REF}
    else
        B2_INDEX_GENO1=${B2_INDEX_REF}
    fi
fi
if [[ ${ID_GENO2} == 'REF' ]]
then
    ID_REF=$(basename ${REF_GENO%.fa*})
    if [[ -z ${B2_INDEX_REF} ]]
    then
        B2_INDEX_GENO2=${INDEXES}/${ID_REF}
    else
        B2_INDEX_GENO2=${B2_INDEX_REF}
    fi
fi

#- Indexes ?
if [[ -z ${B2_INDEX_GENO1} ]]; then B2_INDEX_GENO1=${INDEXES}/${ID_GENO1};fi
if [[ ! -e ${B2_INDEX_GENO1}.rev.2.bt2 ]]
then
    echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_GENO1}. Exit." 1>&2
    echo "$0: ERROR - HINT: Generate genomes and indexes using build_genomes_indexes.sh." 1>&2
    exit 1
fi
if [[ -z ${B2_INDEX_GENO2} ]]; then B2_INDEX_GENO2=${INDEXES}/${ID_GENO2};fi
if [[ ! -e ${B2_INDEX_GENO2}.rev.2.bt2 ]]
then
    echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_GENO2}. Exit." 1>&2
    echo "$0: ERROR - HINT: Generate genomes and indexes using build_genomes_indexes.sh." 1>&2
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
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_GENO1} -U ${FQ_READS_F} 2> ${BAM_OUT}/${ID_OUTBAM}.b2logs | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_GENO1}.bam
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_GENO2} -U ${FQ_READS_F} 2> ${BAM_OUT}/${ID_OUTBAM}.b2logs | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_GENO2}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -p ${BAM_OUT}/${ID_GENO1}.bam -m ${BAM_OUT}/${ID_GENO2}.bam -o ${BAM_OUT} -n ${ID_OUTBAM}
else
    # Paired-end mapping
    # - Mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_GENO1} -1 ${FQ_READS_F} -2 ${FQ_READS_R} 2> ${BAM_OUT}/${ID_OUTBAM}.b2logs | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_GENO1}.bam
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_GENO2} -1 ${FQ_READS_F} -2 ${FQ_READS_R} 2> ${BAM_OUT}/${ID_OUTBAM}.b2logs | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_GENO2}.bam
    # - Selection of best alignments
    ${PYTHON} ${MERGE_ALIGN} -p ${BAM_OUT}/${ID_GENO1}.bam -m ${BAM_OUT}/${ID_GENO2}.bam -s 2 -o ${BAM_OUT} -n ${ID_OUTBAM}
fi

exit 0
