#!/usr/bin/env bash
# Author(s) : Laur (ne Syx
# Contact : laurene.syx@curie.fr
# Comment(s) : From ChIPseq input to filtered VCF file 

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-b"" <ChIPseq input BAM file>"
    echo -e "-v"" <VCF file>"
    echo -e "-o"" <Output directory>"
    echo -e "-h"" <help>"
    exit
}

#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-c) CONFIG=$2; shift;;
    (-b) BAM_FILE=$2; shift;;
    (-v) VCF_FILE=$2; shift;;
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
source /bioinfo/local/curie/ngs-data-management/virtualenvs/hgvs0.3.7/bin/activate

#### Main #### --------------------------------------------------------------------------

# SAMTools mpileup

${SAMTOOLS_DIR}/samtools mpileup -B -Q 0 -d 100000 -l ${FULL_VCF} -f ${REF_GENO} ${BAM_FILE} > ${OUT_DIR}/${OUT_NAME}.pileup

# Count for each base of the VCF file (check_variants.py)

echo -e "${OUT_NAME}\t${OUT_DIR}/${OUT_NAME}.pileup" > ${OUT_DIR}/${OUT_NAME}.conf
${PYTHON} ${SCRIPTS_DIR}/clinTools/check_variants.py ${OUT_DIR}/${OUT_NAME}.conf ${REF_GENO} > ${OUT_DIR}/${OUT_NAME}_covInfo.data

# Allele-specific count for each SNP

awk '(NR>1){OFS="\t";print $2,$3-1,$3,$4,$5,$6,$7,$8,$9,$2,$3-1,$3,"NoName",0,"+"}' ${OUT_DIR}/${OUT_NAME}_covInfo.data > ${OUT_DIR}/${OUT_NAME}_covGene.data
${PYTHON} ${SCRIPTS_DIR}/alleleCount.py --r1 ${OUT_DIR}/${OUT_NAME}_covGene.data --strand 0 -v ${SNP_FILE} -1 ${ID_GENO1} -2 ${ID_GENO2} -o ${OUT_DIR} -n ${OUT_NAME}