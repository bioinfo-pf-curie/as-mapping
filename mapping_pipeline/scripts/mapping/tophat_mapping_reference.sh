#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   Script to map reads on reference genome using Tophat
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
	exit
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

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_reference_tophat
mkdir -p ${BAM_OUT}/

# Name of the reference genome for indexes
ID_REF=$(basename ${REF_GENO%.fa*})
# Name for the output bam file
ID_OUTBAM=${OUT_NAME}_reference

# Checking input parameters for mapping
#  Indexes ?
if [[ -z ${B2_INDEX_REF} ]]; then B2_INDEX_REF=${INDEXES}/${ID_REF};fi
if [[ ! -e ${B2_INDEX_REF}.rev.2.bt2 ]]
then
	echo "$0: ERROR - Missing Bowtie2 indexes for ${ID_REF}. Exit." 1>&2
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
#  SNP file ?
if [[ -z ${SNP_FILE} ]]; then SNP_FILE=${FASTA_OUT}/all_SNPs_${PATERNAL}_${MATERNAL}.txt.gz;fi
if [[ ! -e ${SNP_FILE} ]]
then
    echo "$0: ERROR - Missing SNP file for markAllelicStatus." 1>&2
    exit 1
fi
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
    #${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_REF} ${FQ_READS_F}
    mv ${BAM_OUT}/tophat_out/accepted_hits.bam ${BAM_OUT}/${ID_OUTBAM}.bam
    ${MARKALLELICSTATUS} -r -f -i ${BAM_OUT}/${ID_OUTBAM}.bam -s ${SNP_FILE} -o ${BAM_OUT}/${ID_OUTBAM}_flagged.bam
else
    # Paired-end mapping
    #${TOPHAT_DIR}/tophat -o ${BAM_OUT}/tophat_out ${TOPHAT_OPT} ${TOPHAT_B2_OPT} ${B2_INDEX_REF} ${FQ_READS_F} ${FQ_READS_R}
    mv ${BAM_OUT}/tophat_out/accepted_hits.bam ${BAM_OUT}/${ID_OUTBAM}.bam
    ${MARKALLELICSTATUS} -r -f --paired -i ${BAM_OUT}/${ID_OUTBAM}.bam -s ${SNP_FILE} -o ${BAM_OUT}/${ID_OUTBAM}_flagged.bam
fi


# Filter [TO BE DONE] using markAllelicStatus

exit 0 
