#!/bin/bash

## Nicolas Servant
## Create a masked index for alignment


BEDTOOLS_PATH=/bioinfo/local/build/BEDTools/BEDTools_2.21.0/bin
BOWTIE2_PATH=/bioinfo/local/build/bowtie2/bowtie2-2.2.5/


## Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-r"" <reference name>"
    echo -e "-f"" <per chromosome fasta directory>"
    echo -e "-v"" <VCF/BED file with SNPs>"
    echo -e "-o"" <output directory>"
    echo -e "-h"" <help>"
    exit
}

OUTPUT=$(pwd)
while [ $# -gt 0 ]
do
    case "$1" in
	(-r) ORG=$2; shift;;
	(-f) FASTA=$2; shift;;
	(-v) VCF=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

#1- create a masked fasta from VCF file
if [[ ! -d ${FASTA} || ! -e ${VCF} ]]; then echo "Error: Input files not found. Exit"; exit; fi
if [[ -z ${ORG} ]]; then echo "Error : Please specify a reference name. Exit"; exit; fi

echo "Masking fasta files ..."
mkdir -p ${OUTPUT}maskfasta
for fa in `ls $FASTA*.fa`; do 
    if [[ ! $fa =~ random && ! $fa =~ chrMT ]]; then
	echo ${fa}
	OUT=`basename $fa | sed -e 's/.fa/_masked.fa/'`
	${BEDTOOLS_PATH}/bedtools maskfasta -fi ${fa} -bed ${VCF} -fo ${OUTPUT}maskfasta/${OUT}
    fi
done

cat ${OUTPUT}maskfasta/*.fa >> ${OUTPUT}${ORG}.fa
rm -r ${OUTPUT}maskfasta/

#2- create the index from the masked fasta file 
echo "Generating bowtie2 indexes ..."
mkdir -p ${OUTPUT}bowtie2_indexes
${BOWTIE2_PATH}bowtie2-build -f ${OUTPUT}${ORG}.fa ${OUTPUT}bowtie2_indexes/${ORG}
