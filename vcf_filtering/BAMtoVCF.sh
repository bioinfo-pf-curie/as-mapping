#!/usr/bin/env bash
# Author(s) : Laurene Syx
# Contact : laurene.syx@curie.fr
# Comment(s) : From ChIPseq input to filtered VCF file 

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-b"" <ChIPseq input BAM file>"
    echo -e "-v"" <SNP list file>"
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
    (-v) SNP_FILE=$2; shift;;
    (-o) OUT_DIR=$2; shift;;
    (-n) OUT_NAME=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: ERROR - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

source /bioinfo/local/curie/ngs-data-management/virtualenvs/hgvs0.3.7/bin/activate

#### Paths #### --------------------------------------------------------------------------

BEDTOOLS_DIR=/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin
SAMTOOLS_DIR=/bioinfo/local/build/samtools/samtools-1.1/bin
PYTHON=/bioinfo/local/build/python/python-2.7.9/bin/python2.7
SCRIPTS_DIR=/bioinfo/users/lsyx/GitLab/as-mapping/vcf_filtering/scripts

ID_GENO1=   #REF
ID_GENO2=   #ALT

#### Main #### --------------------------------------------------------------------------

mkdir -p ${OUT_DIR}

# SNP file to BED / reformat SNP file
zcat ${SNP_FILE} | awk 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT"}{OFS="\t"; split($5,snp,"/");print $2,$3,snp[1],snp[2]}' > ${OUT_DIR}/${OUT_NAME}_SNPlist.txt
awk -F '\t' '($1!~/^#/){OFS="\t"; print $1,$2}' ${OUT_DIR}/${OUT_NAME}_SNPlist.txt > ${OUT_DIR}/${OUT_NAME}_SNPlist.pos

# BAM sort 
${SAMTOOLS_DIR}/samtools sort -T ${OUT_DIR}/${OUT_NAME}_sorted_tmp -o ${OUT_DIR}/${OUT_NAME}_sorted.bam ${BAM_FILE}

# SAMTools mpileup
${SAMTOOLS_DIR}/samtools mpileup -B -Q 0 -d 100000 -l ${OUT_DIR}/${OUT_NAME}_SNPlist.pos -f ${REF_GENO} ${OUT_DIR}/${OUT_NAME}_sorted.bam > ${OUT_DIR}/${OUT_NAME}.pileup

# Count for each base of the VCF file (check_variants.py)
echo -e "${OUT_NAME}\t${OUT_DIR}/${OUT_NAME}.pileup" > ${OUT_DIR}/${OUT_NAME}.conf
${PYTHON} ${SCRIPTS_DIR}/clinTools/check_variants.py ${OUT_DIR}/${OUT_NAME}.conf ${REF_GENO} > ${OUT_DIR}/${OUT_NAME}_covInfo.data

# Allele-specific count for each SNP
awk -F '\t' '(NR>1){OFS="\t";print $2,$3-1,$3,$6,$7,$8,$9}' ${OUT_DIR}/${OUT_NAME}_covInfo.data > ${OUT_DIR}/${OUT_NAME}_covInfo_tmp.data
mv ${OUT_DIR}/${OUT_NAME}_covInfo_tmp.data ${OUT_DIR}/${OUT_NAME}_covInfo.data
${PYTHON} ${SCRIPTS_DIR}/alleleCount.py -f ${OUT_DIR}/${OUT_NAME}_covInfo.data -v ${OUT_DIR}/${OUT_NAME}_SNPlist.txt -1 ${ID_GENO1} -2 ${ID_GENO2} -o ${OUT_DIR}/ -n ${OUT_NAME}_SNPInfo

# Heterozygous SNPs
awk -F '\t' '(NR>1){OFS="\t"; totalCount=$4+$9; if(totalCount>=10){print $1, $2, $3, totalCount, $6/totalCount, $8/totalCount, $9/totalCount}}' ${OUT_DIR}/${OUT_NAME}_SNPInfo.data  | awk -F '\t' '($1!="chrY" && $1!="chrM" && $7<=0.1 && $5>=0.25 && $5<=0.75 && $6>=0.25 && $6<=0.75){OFS="\t"; print $1,$2,$3}' > ${OUT_DIR}/${OUT_NAME}.bed

# Generate the clean SNP file
 awk -F '\t' '(NR>1){OFS="\t"; print $1, $2-1, $2, $3, $4}' ${OUT_DIR}/${OUT_NAME}_SNPlist.txt | ${BEDTOOLS_DIR}/intersectBed -a stdin -b ${OUT_DIR}/${OUT_NAME}.bed -wa | awk -F '\t' '{OFS="\t"; sub("chr","",$1); print 1,$1,$3,1,$4"/"$5}' > ${OUT_DIR}/${OUT_NAME}.txt

# Remove all tmp files
rm ${OUT_DIR}/${OUT_NAME}_SNPlist.txt ${OUT_DIR}/${OUT_NAME}_SNPlist.pos ${OUT_DIR}/${OUT_NAME}_sorted.bam ${OUT_DIR}/${OUT_NAME}.pileup ${OUT_DIR}/${OUT_NAME}_covInfo.data ${OUT_DIR}/${OUT_NAME}.conf ${OUT_DIR}/${OUT_NAME}_SNPInfo.data ${OUT_DIR}/${OUT_NAME}.bed 
