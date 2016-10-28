#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Script to reconstruct N-masked genome from two parents in fasta format from VCF files containing SNPs

#### Paths #### ----------------------------------------------------------------------

SNPSPLIT_GEN=/data/annotations/Mouse/personalized_genomes/SNPsplit/SNPsplit_genome_preparation

#### Functions #### ------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-p/--paternal"" <Name of the paternal strain>" 
    echo -e "-m/--maternal"" <Name of the maternal strain>" 
    echo -e "-v/--vcf"" <VCF file>" 
    echo -e "-r/--reference"" <Directory with reference chromosomes (Ensembl format)>" 
    echo -e "-o/--outdir"" <Output directory>" 
    echo -e "[-g/--snpsplit_gen]"" <SNPsplit_genome_preparation path>" 
    echo -e "-h/--help"" <Help>"
    exit 0
}

#### Read arguments #### -------------------------------------------------------------

ARG_LIST=($@)

# Long arguments parsing
for arg in ${ARG_LIST[@]}
do
    delim=""
    case "${arg}" in
        --help) args="${args}-h ";;
        --paternal) args="${args}-p ";;
        --maternal) args="${args}-m ";;
        --vcf) args="${args}-v ";;
        --reference) args="${args}-r ";;
        --outdir) args="${args}-o ";;
        --snpsplit_gen) args="${args}-g ";;
        *) [[ "${arg:0:1}" == "-" ]] || delim="\""
            args="${args}${delim}${arg}${delim} ";;
    esac
done

eval set -- ${args}

# Short arguments parsing
while getopts :p:m:v:r:o:g:h option
do
    case "${option}" in
        p) PATERNAL=${OPTARG};;
        m) MATERNAL=${OPTARG};;
        v) VCF=${OPTARG};;
        r) REF_DIR=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        g) SNPSPLIT_GEN=${OPTARG};;
        h) usage;;
        \?) echo "$0: ERROR - Invalid option: - ${OPTARG}" 1>&2; exit 1;;
        :) echo "$0: ERROR - Option -${OPTARG} requires an argument." 1>&2; exit 1;;
    esac
done

#### Main #### -----------------------------------------------------------------------

mkdir -p ${OUT_DIR}

# -- Run SNPsplit_genome_preparation to generate parental and N-masked genomes
echo "$0: Generating N-masked genomes of $PATERNAL and $MATERNAL ..."

# In case one of the genotype is from reference genome (REF)
if [[ ${PATERNAL} == 'REF' ]]
then
    STRAINS="--strain ${MATERNAL}"
    PREFIX=${MATERNAL}
elif [[ ${MATERNAL} == 'REF' ]]
then
    STRAINS="--strain ${PATERNAL}"
    PREFIX=${PATERNAL}
else # Hybrid genome with strains different from reference genome
    STRAINS="--strain ${PATERNAL} --strain2 ${MATERNAL}"
    PREFIX=${PATERNAL}_${MATERNAL}
fi

${SNPSPLIT_GEN} --nmasking ${STRAINS} --reference_genome ${REF_DIR} --vcf_file ${VCF}

# -- Concatenation of the N-masked genome
for i in `ls -d --color=never ${PREFIX}*_N-masked/*`
do
    if [[ $i =~ MT ]]
    then
        sed 's/>MT/>chrM/' ${i} >> ${OUT_DIR}/N-masked_${PATERNAL}_${MATERNAL}.fa
    else
	    sed 's/>/>chr/' ${i} >> ${OUT_DIR}/N-masked_${PATERNAL}_${MATERNAL}.fa
    fi
done

# -- Save and delete files

# Reports
mkdir -p ${OUT_DIR}/SNPsplit_reports
mv ${PREFIX}_*report.txt ${OUT_DIR}/SNPsplit_reports

# SNPs and delete specific files
if [[ ${PATERNAL} == 'REF' ]] || [[ ${MATERNAL} == 'REF' ]] 
then
    gunzip -c all_SNPs_${PREFIX}_*.txt.gz | awk '{if ($2=="MT") $2="chrM";else $2="chr"$2; print}' OFS='\t' > ${OUT_DIR}/all_SNPs_${PATERNAL}_${MATERNAL}.txt
    rm all_SNPs_${PREFIX}_*.txt.gz
    rm -r SNPs_${PREFIX} ${PREFIX}_N-masked
else
    awk '{if ($2=="MT") $2="chrM";else $2="chr"$2; print}' OFS='\t' all_${MATERNAL}_SNPs_${PATERNAL}_*.txt > ${OUT_DIR}/all_SNPs_${PATERNAL}_${MATERNAL}.txt
    # Delete specific files
    rm all_${MATERNAL}_SNPs_${PATERNAL}_*.txt
    rm -r ${PATERNAL}_${MATERNAL}_*_full_sequence ${PATERNAL}_${MATERNAL}_*_N-masked ${PATERNAL}_${MATERNAL}_SNPs_in_common.*.txt
    for GENO in $PATERNAL $MATERNAL
    do
        rm -r ${GENO}_full_sequence ${GENO}_N-masked SNPs_${GENO}
        rm *_${GENO}_*.gz ${GENO}_*report.txt ${GENO}_specific_*.txt
    done
fi

# Compress SNP file
gzip ${OUT_DIR}/all_SNPs_${PATERNAL}_${MATERNAL}.txt

exit 0
