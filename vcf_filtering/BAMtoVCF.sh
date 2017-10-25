#!/usr/bin/env bash
# Author(s) : Laurene Syx
# Contact : laurene.syx@curie.fr
# Comment(s) : From ChIPseq input to filtered VCF file 

set -o pipefail
set -o errexit

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-b"" <ChIPseq input BAM file>"
    echo -e "-v"" <SNP list file>"
    echo -e "-f"" <Reference Fasta file>"
    echo -e "-o"" <Output directory>"
    echo -e "-h"" <help>"
    exit
}

#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-b) BAM_FILE=$2; shift;;
    (-v) SNP_FILE=$2; shift;;
    (-f) REF_GENO=$2; shift;;
    (-o) OUT_DIR=$2; shift;;
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

CUR_DIR=$(realpath $(dirname $0))
prefix=$(basename $BAM_FILE | sed -e 's/.bam//')

ID_GENO1=ref   #REF
ID_GENO2=alt   #ALT

## OPTIONS
min_cov=5
ht_min=0.2
ht_max=0.8
max_conta=0.1

#### Main #### --------------------------------------------------------------------------

mkdir -p ${OUT_DIR}

# SNP file to BED / reformat SNP file
echo -e "Prepare SNPs file ..."
zcat ${SNP_FILE} | awk 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT"}{OFS="\t"; split($5,snp,"/");print "chr"$2,$3,snp[1],snp[2]}' > ${OUT_DIR}/${prefix}_SNPlist.txt
awk -F '\t' '($1!~/^#/){OFS="\t"; print $1,$2}' ${OUT_DIR}/${prefix}_SNPlist.txt > ${OUT_DIR}/${prefix}_SNPlist.pos

# BAM sort 
is_sorted=$( ${SAMTOOLS_DIR}/samtools view -H ${BAM_FILE} | awk 'BEGIN{r=0}$0~"SO:unsorted"{r=1}END{print r}' )
if [ $is_sorted == 1 ]; then
    echo -e "Sort BAM file ..."
    ${SAMTOOLS_DIR}/samtools sort -@ 4 -T ${OUT_DIR}/${prefix}_sorted_tmp -o ${OUT_DIR}/${prefix}_sorted.bam ${BAM_FILE}
else
    echo -e "BAM file appears to be already sorted by coordinate ..."
    ln -f -s ${BAM_FILE} ${OUT_DIR}/${prefix}_sorted.bam 
fi

# SAMTools mpileup
echo -e "Run samtools mpileup ..."
${SAMTOOLS_DIR}/samtools mpileup -B -Q 0 -d 100000 -l ${OUT_DIR}/${prefix}_SNPlist.pos -f ${REF_GENO} ${OUT_DIR}/${prefix}_sorted.bam > ${OUT_DIR}/${prefix}.pileup

# Count for each base of the VCF file (check_variants.py)
echo -e "Check variants ..."
echo -e "${prefix}\t${OUT_DIR}/${prefix}.pileup" > ${OUT_DIR}/${prefix}.conf
${PYTHON} ${CUR_DIR}/scripts/clinTools/check_variants.py ${OUT_DIR}/${prefix}.conf ${REF_GENO} > ${OUT_DIR}/${prefix}_covInfo.data

# Allele-specific count for each SNP
echo -e "Allele-specific counts and filtering ..."
awk -F '\t' '(NR>1){OFS="\t";print $2,$3-1,$3,$6,$7,$8,$9}' ${OUT_DIR}/${prefix}_covInfo.data > ${OUT_DIR}/${prefix}_covInfo_tmp.data
mv ${OUT_DIR}/${prefix}_covInfo_tmp.data ${OUT_DIR}/${prefix}_covInfo.data
${PYTHON} ${CUR_DIR}/scripts/alleleCount.py -f ${OUT_DIR}/${prefix}_covInfo.data -v ${OUT_DIR}/${prefix}_SNPlist.txt -1 ${ID_GENO1} -2 ${ID_GENO2} -o ${OUT_DIR}/ -n ${prefix}_SNPInfo

# Heterozygous SNPs - 0.25/0.75 for expected alleles + others < 10% + total >= 10
awk -F '\t' -v mincounts=$min_cov '(NR>1){OFS="\t"; totalCount=$4+$9; if(totalCount>=mincounts && $4 > 0){print $1, $2, $3, totalCount, $6/$4, $8/$4, $9/totalCount}}' ${OUT_DIR}/${prefix}_SNPInfo.data  | \
    awk -F '\t' -v htmin=$ht_min -v htmax=$ht_max -v maxconta=$max_conta '($1!="chrY" && $1!="chrM" && $7<=maxconta && $5>=htmin && $5<=htmax){OFS="\t"; print $1,$2,$3}' > ${OUT_DIR}/${prefix}.bed

# Generate the clean SNP file
 awk -F '\t' '(NR>1){OFS="\t"; print $1, $2-1, $2, $3, $4}' ${OUT_DIR}/${prefix}_SNPlist.txt | ${BEDTOOLS_DIR}/intersectBed -a stdin -b ${OUT_DIR}/${prefix}.bed -wa | \
     awk -F '\t' '{OFS="\t"; sub("chr","",$1); print 1,$1,$3,1,$4"/"$5}' > ${OUT_DIR}/${prefix}.txt

# Remove all tmp files
##rm ${OUT_DIR}/${prefix}_SNPlist.txt ${OUT_DIR}/${prefix}_SNPlist.pos ${OUT_DIR}/${prefix}_sorted.bam ${OUT_DIR}/${prefix}.pileup ${OUT_DIR}/${prefix}_covInfo.data ${OUT_DIR}/${prefix}.conf ${OUT_DIR}/${prefix}_SNPInfo.data ${OUT_DIR}/${prefix}.bed 
