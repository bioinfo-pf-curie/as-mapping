#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   Script to map reads on reference genome using Bowtie2
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

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_reference
mkdir -p ${BAM_OUT}/

# Get the name of the reference genome to use already existing indexes (might change this part to work directly on my own indexes)
ID_REF=$(basename ${REF_GENO})
ID_REF=${ID_REF%.fa}

# -- Checking input parameters for mapping
#  Indexes ?
if [[ -z ${B2_INDEX_REF} ]]; then B2_INDEX_REF=${INDEXES}/${ID_REF}
if [[ ! -e ${B2_INDEX_REF}.rev.2.bt2 ]]
then
	echo "$0: ERROR - Missing Bowtie2 indexes in ${INDEXES}. Exit." 1>&2
    exit 1
fi

#  Reads ?
if [[ ! -e ${FQ_READS_F} ]]
then
	echo "$0: ERROR - Missing FASTQ reads. Exit." 1>&2
    exit 1
fi
SINGLE_END = true
if [[ -e ${FQ_READS_R} ]]; then SINGLE_END = false;fi


# -- Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_REF} -U ${FQ_READS_F} | ${SAMTOOLS} view -bS - > ${BAM_OUT}${ID_REF}.bam
else
    # Paired-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_REF} -1 ${FQ_READS_F} -2 ${FQ_READS_R} | ${SAMTOOLS} view -bS - > ${BAM_OUT}${ID_REF}.bam
fi


# -- Filter [TO BE DONE]

exit 0
