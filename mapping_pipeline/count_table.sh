#!/usr/bin/env bash
# Author(s) : Laurène Syx
# Contact : laurene.syx@curie.fr
# Comment(s) : Create the count table

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-b"" <Flagged BAM file>"
    echo -e "-o"" <Output directory for the count files>"
    echo -e "-n"" <Output name for the count files>"
    echo -e "-h"" <help>"
    exit
}


#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-c) CONFIG=$2; shift;;
    (-b) BAM_FILE=$2; shift;;
    (-o) OUT_DIR=$2; shift;;
    (-n) OUT_NAME=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: ERROR - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

if [[ -z ${CONFIG} ]]; then echo "$0: ERROR - you need to specify a CONFIG file (-c). Exit." 1>&2; exit 1; fi
source ${CONFIG}
source ${PIPELINE_PATH}/includes/path_fct.inc

#### Main #### --------------------------------------------------------------------------

# Create output directory
mkdir -p ${OUT_DIR}/count

# Separation in 2 parental specific files: Parent1 and Parent2
${SAMTOOLS_DIR}/samtools view -h ${BAM_FILE} | grep -E '^@|XX:Z:G1' | ${SAMTOOLS_DIR}/samtools view -b - > ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO1}.bam
${SAMTOOLS_DIR}/samtools view -h ${BAM_FILE} | grep -E '^@|XX:Z:G2' | ${SAMTOOLS_DIR}/samtools view -b - > ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO2}.bam

# Count for each gene
${FEATURECOUNTS_DIR} ${FEATURECOUNTS_OPT} -a ${REF_GTF} -o ${OUT_DIR}/count/${OUT_NAME}_count.txt ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO1}.bam ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO2}.bam ${BAM_FILE}
awk -F '\t' -v parent1=${ID_GENO1} -v parent2=${ID_GENO2} -v threshold=${READS_THRESHOLD} 'BEGIN{OFS="\t"; print "Gene", parent1, parent2, "allelicRatio_"parent1"/all", "allReads"}{if($1!~/^#/ && $1!="Geneid"){if($7+$8>=threshold){as=$7/($7+$8)}else{as="NA"};print $1,$7,$8,as,$9}}' ${OUT_DIR}/count/${OUT_NAME}_count.txt >  ${OUT_DIR}/count/${OUT_NAME}_allelicRatio.txt


# Delete all temporary BAM files
/bin/rm ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO1}.bam ${OUT_DIR}/count/${OUT_NAME}_${ID_GENO2}.bam 

exit 0 
