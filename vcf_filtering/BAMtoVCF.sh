#!/usr/bin/env bash
# Author(s) : Laurene Syx
# Contact : laurene.syx@curie.fr
# Modified : nicolas.servant@curie.fr - 27-10-17
# Comment(s) : From ChIPseq input to filtered VCF file 

##
## BAMtoVCF.sh
## Combine BAM file with variant information to select SNPs that can be used for allele-specific mapping
##
## bash BAMtoVCF.sh -i $BAM -s $SNP -f /data/annotations/Mouse/mm9/mm9.fa -o ./test/ -d 5 -t 0.25 -c 0.2 -a -v


set -o pipefail
set -o errexit

#### Function #### ----------------------------------------------------------------------
function usage {
    echo -e "usage : $0 -i BAM_FILE -s SNP_FILE -f REF_GENO -o OUTPUT_DIR [-d MIN_COV] [-t HT_TH] [-c MAX_CONTA] [-a] [-h] [-v]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo
    echo -e "Help : $0"
    echo -e "-i"" <ChIPseq input BAM file>"
    echo -e "-s"" <SNP list file (zipped) - see SNPsplit format ; ID CHR POS 1 REF/ALT>"
    echo -e "-f"" <Reference Fasta file>"
    echo -e "-o"" [<Output directory. Default is current folder>]"
    echo -e "-a"" [<Include abnormal pairs in the mpileup. Default=No>]"
    echo -e "-d"" [<Minimum read depth to consider a SNP. Default=10>]"
    echo -e "-t"" [<Allelic ratio range (0.5 +- threshold) to consider a SNP as heterozygote. Default=0.25>]"
    echo -e "-c"" [<Maximum fraction of unexpected bases per SNP. Default=0.1>]"
    echo -e "-v"" <verbose>"
    echo -e "-h"" <help>"
    exit
}

#### Parameters #### --------------------------------------------------------------------
KEEP_ABNORMAL_PAIRS=0
VERBOSE=0
MIN_COV=10
HT_TH=0.25
MAX_CONTA=0.1
ODIR="./"

if [ $# -lt 1 ];
then
    usage
    exit
fi

while getopts ":i:s:f:o:d:t:c:avh" OPT
do
    case $OPT in
	i) BAM_FILE=$OPTARG;;
        s) SNP_FILE=$OPTARG;;
        f) REF_GENO=$OPTARG;;
        o) ODIR=$OPTARG;;
        a) KEEP_ABNORMAL_PAIRS=1;;
        d) MIN_COV=$OPTARG;;
        t) HT_TH=$OPTARG;;
        c) MAX_CONTA=$OPTARG;;
        v) VERBOSE=1 ;;
        h) help ;;
        \?)
             echo "Invalid option: -$OPTARG" >&2
             usage
             exit 1
             ;;
         :)
             echo "Option -$OPTARG requires an argument." >&2
             usage
             exit 1
             ;;
    esac
done

source /bioinfo/local/curie/ngs-data-management/virtualenvs/hgvs0.3.7/bin/activate

if [[ ! -e $BAM_FILE || ! -e $SNP_FILE ]]
then
    >&2 echo -e "Error - Input file(s) not found !"
    exit 1
fi

#### Paths #### --------------------------------------------------------------------------

BEDTOOLS_DIR=/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin
SAMTOOLS_DIR=/bioinfo/local/build/samtools/samtools-1.1/bin
PYTHON=/bioinfo/local/build/python/python-2.7.9/bin/python2.7

CUR_DIR=$(realpath $(dirname $0))
prefix=$(basename $BAM_FILE | sed -e 's/.bam//')

id_geno1=ref   #REF
id_geno2=alt   #ALT
ht_min=$(echo "0.5 - $HT_TH" | bc)
ht_max=$(echo "0.5 + $HT_TH" | bc)

#### Main #### --------------------------------------------------------------------------
if [ $VERBOSE == 1 ]; then
    echo -e "## BAM2VCF"
    echo -e "## INPUT=$BAM_FILE"
    echo -e "## SNP_FILE=$SNP_FILE"
    echo -e "## REFERENCE=$REF_GENO"
    echo -e "## OUTPUT_DIRECTORY=$ODIR"
    echo -e "## "
    echo -e "## KEEP_ABNORMAL_PAIRS=$KEEP_ABNORMAL_PAIRS"
    echo -e "## MIN_DEPTH=$MIN_COV"
    echo -e "## HETERO_RANGE=[$ht_min - $ht_max]"
    echo -e "## MAX_CONTA=$MAX_CONTA"
    echo -e ""
fi

mkdir -p ${ODIR}

# SNP file to BED / reformat SNP file
if [ $VERBOSE == 1 ]; then echo -e "## Prepare SNPs file ..."; fi
zcat ${SNP_FILE} | awk 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT"}{OFS="\t"; split($5,snp,"/");print "chr"$2,$3,snp[1],snp[2]}' > ${ODIR}/${prefix}_SNPlist.txt
awk -F '\t' '($1!~/^#/){OFS="\t"; print $1,$2}' ${ODIR}/${prefix}_SNPlist.txt > ${ODIR}/${prefix}_SNPlist.pos

# BAM sort 
is_sorted=$( ${SAMTOOLS_DIR}/samtools view -H ${BAM_FILE} | awk 'BEGIN{r=0}$0~"SO:unsorted"{r=1}END{print r}' )
if [ $is_sorted == 1 ]; then
    if [ $VERBOSE == 1 ]; then echo -e "## Sort BAM file ..."; fi
    ${SAMTOOLS_DIR}/samtools sort -@ 4 -T ${ODIR}/${prefix}_sorted_tmp -o ${ODIR}/${prefix}_sorted.bam ${BAM_FILE}
else
    if [ $VERBOSE == 1 ]; then echo -e "## BAM file appears to be already sorted by coordinate ..."; fi
    ln -f -s ${BAM_FILE} ${ODIR}/${prefix}_sorted.bam 
fi

# SAMTools mpileup
if [ $KEEP_ABNORMAL_PAIRS == 1 ]; then
    if [ $VERBOSE == 1 ]; then echo -e "## Run samtools mpileup with abnormal pairs ..."; fi
    ${SAMTOOLS_DIR}/samtools mpileup -A -B -Q 0 -d 100000 -l ${ODIR}/${prefix}_SNPlist.pos -f ${REF_GENO} ${ODIR}/${prefix}_sorted.bam > ${ODIR}/${prefix}.pileup
else
    if [ $VERBOSE == 1 ]; then echo -e "## Run samtools mpileup ..."; fi
    ${SAMTOOLS_DIR}/samtools mpileup -B -Q 0 -d 100000 -l ${ODIR}/${prefix}_SNPlist.pos -f ${REF_GENO} ${ODIR}/${prefix}_sorted.bam > ${ODIR}/${prefix}.pileup
fi

# Count for each base of the VCF file (check_variants.py)
if [ $VERBOSE == 1 ]; then echo -e "## Check variants ..."; fi
echo -e "${prefix}\t${ODIR}/${prefix}.pileup" > ${ODIR}/${prefix}.conf
${PYTHON} ${CUR_DIR}/scripts/clinTools/check_variants.py ${ODIR}/${prefix}.conf ${REF_GENO} > ${ODIR}/${prefix}_covInfo.data

# Allele-specific count for each SNP
if [ $VERBOSE == 1 ]; then echo -e "## Allele-specific counts and filtering ..."; fi
awk -F '\t' '(NR>1){OFS="\t";print $2,$3-1,$3,$6,$7,$8,$9}' ${ODIR}/${prefix}_covInfo.data > ${ODIR}/${prefix}_covInfo_tmp.data
mv ${ODIR}/${prefix}_covInfo_tmp.data ${ODIR}/${prefix}_covInfo.data
${PYTHON} ${CUR_DIR}/scripts/alleleCount.py -f ${ODIR}/${prefix}_covInfo.data -v ${ODIR}/${prefix}_SNPlist.txt -1 ${id_geno1} -2 ${id_geno2} -o ${ODIR}/ -n ${prefix}_SNPInfo

# Heterozygous SNPs - 0.25/0.75 for expected alleles + others < 10% + total >= 10
awk -F '\t' -v mincounts=$MIN_COV '(NR>1){OFS="\t"; totalCount=$4+$9; if(totalCount>=mincounts && $4 > 0){print $1, $2, $3, totalCount, $6/$4, $8/$4, $9/totalCount}}' ${ODIR}/${prefix}_SNPInfo.data  | \
    awk -F '\t' -v htmin=$ht_min -v htmax=$ht_max -v maxconta=$MAX_CONTA '($1!="chrY" && $1!="chrM" && $7<maxconta && $5>=htmin && $5<=htmax){OFS="\t"; print $1,$2,$3}' > ${ODIR}/${prefix}.bed

# Generate the clean SNP file
 awk -F '\t' '(NR>1){OFS="\t"; print $1, $2-1, $2, $3, $4}' ${ODIR}/${prefix}_SNPlist.txt | ${BEDTOOLS_DIR}/intersectBed -a stdin -b ${ODIR}/${prefix}.bed -wa | \
     awk -F '\t' '{OFS="\t"; sub("chr","",$1); print 1,$1,$3,1,$4"/"$5}' > ${ODIR}/${prefix}.txt

# Remove all tmp files
rm ${ODIR}/${prefix}_SNPlist.txt ${ODIR}/${prefix}_SNPlist.pos ${ODIR}/${prefix}_sorted.bam ${ODIR}/${prefix}.pileup ${ODIR}/${prefix}_covInfo.data ${ODIR}/${prefix}.conf ${ODIR}/${prefix}_SNPInfo.data ${ODIR}/${prefix}.bed 

if [ $VERBOSE == 1 ]; then echo -e "## done !" ; fi