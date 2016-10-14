#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
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
BAM_OUT=${OUT_DIR}/mapping_reference_bowtie2
mkdir -p ${BAM_OUT}/

# Name of the reference genome for indexes
ID_REF=$(basename ${REF_GENO%.fa*})
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_reference

# -- Checking input parameters for mapping
#  Indexes ?
if [[ -z ${B2_INDEX_REF} ]]; then B2_INDEX_REF=${INDEXES}/${ID_REF};fi
if [[ ! -e ${B2_INDEX_REF}.rev.2.bt2 ]]
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
#  SNP file ?
if [[ -z ${SNP_FILE} ]]; then SNP_FILE=${FASTA_OUT}/all_SNPs_${PATERNAL}_${MATERNAL}.txt.gz;fi
if [[ ! -e ${SNP_FILE} ]]
then
    echo "$0: ERROR - Missing SNP file for markAllelicStatus." 1>&2
    exit 1
fi

# -- Mapping
if [[ $SINGLE_END == true ]]
then
    # Single-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_REF} -U ${FQ_READS_F} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    ${MARKALLELICSTATUS} -r -i ${BAM_OUT}/${ID_OUTBAM}.bam -s ${SNP_FILE} -o ${BAM_OUT}/${ID_OUTBAM}_flagged.bam
else
    # Paired-end mapping
    ${BOWTIE2_DIR}/bowtie2 ${B2_OPT} -x ${B2_INDEX_REF} -1 ${FQ_READS_F} -2 ${FQ_READS_R} | ${SAMTOOLS} view -bS - > ${BAM_OUT}/${ID_OUTBAM}.bam
    ${MARKALLELICSTATUS} -r --paired -i ${BAM_OUT}/${ID_OUTBAM}.bam -s ${SNP_FILE} -o ${BAM_OUT}/${ID_OUTBAM}_flagged.bam
fi

exit 0
