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
    echo -e "[-a/--parental]"" <Also build parental genomes>" 
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
        --parental) args="${args}-a ";;
        *) [[ "${arg:0:1}" == "-" ]] || delim="\""
            args="${args}${delim}${arg}${delim} ";;
    esac
done

eval set -- ${args}

# Short arguments parsing
while getopts :p:m:v:r:o:g:ah option
do
    case "${option}" in
        p) PATERNAL=${OPTARG};;
        m) MATERNAL=${OPTARG};;
        v) VCF=${OPTARG};;
        r) REF_DIR=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        g) SNPSPLIT_GEN=${OPTARG};;
        a) BUILD_PAR=true;;
        h) usage;;
        \?) echo "$0: ERROR - Invalid option: - ${OPTARG}" 1>&2; exit 1;;
        :) echo "$0: ERROR - Option -${OPTARG} requires an argument." 1>&2; exit 1;;
    esac
done

#### Main #### -----------------------------------------------------------------------

mkdir -p ${OUT_DIR}

# -- Run SNPsplit_genome_preparation to generate parental and N-masked genomes
echo "Generating N-masked genomes of $PATERNAL and $MATERNAL ..."

${SNPSPLIT_GEN} --nmasking --strain ${PATERNAL} --strain2 ${MATERNAL} --reference_genome ${REF_DIR} --vcf_file ${VCF}

# -- Concatenation of the N-masked genome
for i in `seq 1 19` X Y MT
do
    if [[ $i == MT ]]
    then
        sed 's/>MT/>chrM/' ${PATERNAL}_${MATERNAL}_*_N-masked/chr${i}.N-masked.fa >> ${OUT_DIR}/N-masked_${PATERNAL}_${MATERNAL}.fa
    else
	    sed 's/>/>chr/' ${PATERNAL}_${MATERNAL}_*_N-masked/chr${i}.N-masked.fa >> ${OUT_DIR}/N-masked_${PATERNAL}_${MATERNAL}.fa
    fi
done

# -- Save and delete files
# Reports
mkdir -p ${OUT_DIR}/SNPsplit_reports
mv ${PATERNAL}_${MATERNAL}_*report.txt ${OUT_DIR}/SNPsplit_reports
# SNPs
gzip all_${MATERNAL}_SNPs_${PATERNAL}_*.txt
mv all_${MATERNAL}_SNPs_${PATERNAL}_*.txt.gz ${OUT_DIR}
# Delete files
rm ${PATERNAL}_${MATERNAL}_SNPs_*.txt
rm -r ${PATERNAL}_${MATERNAL}_*_full_sequence ${PATERNAL}_${MATERNAL}_*_N-masked

# -- [OPTION] Concatenetion of the parental genome in fasta files
if [[ $BUILD_PAR ]]
then
    for GENO in $PATERNAL $MATERNAL
    do
	    for i in `seq 1 19` X Y MT 
	    do
            if [[ $i == MT ]]
            then
                 sed 's/>MT/>chrM/' ${GENO}_full_sequence/chr${i}.SNPs_introduced.fa >> ${OUT_DIR}/${GENO}.fa
		    else
                 sed 's/>/>chr/' ${GENO}_full_sequence/chr${i}.SNPs_introduced.fa >> ${OUT_DIR}/${GENO}.fa
	        fi
        done
        
        # Save and delete files
        # Reports
        mkdir -p ${OUT_DIR}/SNPsplit_reports
        mv ${GENO}_*report.txt ${OUT_DIR}/SNPsplit_reports
        # SNPs
        mv *_${GENO}_*.gz ${OUT_DIR}
        # Delete files
        rm ${GENO}_specific_*.txt
        rm -r SNPs_${GENO} ${GENO}_full_sequence ${GENO}_N-masked
    done    
else
    for GENO in $PATERNAL $MATERNAL
    do
        # Delete files
        rm -r ${GENO}_full_sequence ${GENO}_N-masked SNPs_${GENO}
        rm *_${GENO}_*.gz ${GENO}_*report.txt ${GENO}_specific_*.txt
    done
fi

exit 0
