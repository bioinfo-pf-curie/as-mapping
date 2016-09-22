#!/usr/bin/env bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Script to reconstruct a diploid genome in fasta format from VCF file
#              containing SNPs. Also reconstruct both parental genomes

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
    exit 1
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
        \?) echo "Invalid option: - ${OPTARG}"; exit 1;;
        :) echo "Option -${OPTARG} requires an argument." >&2; exit 1;;
    esac
done


#### Main #### -----------------------------------------------------------------------

mkdir -p $OUT_DIR

# 1- Run SNPsplit_genome_preparation to generate parental genomes
# 2- Concatenation of the parental chromosomes in fasta files
#     - In corresponding parental genomes
#     - In the diploid genome

FASTA_DIPLOID=${PATERNAL}_${MATERNAL}.fa && rm -f ${OUT_DIR}/${FASTA_DIPLOID}

for STRAIN in $PATERNAL $MATERNAL
do
    # -- SNPsplit
    echo "Generating parental genome of ${STRAIN} ..."
    echo "in $OUT_DIR"
    ${SNPSPLIT_GEN} --strain ${STRAIN} --reference_genome ${REF_DIR} --vcf_file ${VCF} --no_nmasking
        
    # -- Concatenation in parental genome
    for i in `seq 1 19` X Y MT
    do
        if [[ $i == MT ]]
        then
            sed 's/>MT/>chrM/' ${STRAIN}_full_sequence/chr${i}.SNPs_introduced.fa >> ${OUT_DIR}/${STRAIN}.fa
        else
            sed 's/>/>chr/' ${STRAIN}_full_sequence/chr${i}.SNPs_introduced.fa >> ${OUT_DIR}/${STRAIN}.fa
        fi
    done

    # -- Add to the diploid genome
    awk -v STRAIN=${STRAIN}  '{if ($1 ~ /^>/) print $1"_"STRAIN; else print}' ${OUT_DIR}/${STRAIN}.fa >> ${OUT_DIR}/${FASTA_DIPLOID}

    # Save and delete files
    # Reports
    mkdir -p ${OUT_DIR}/SNPsplit_reports
    mv *report.txt ${OUT_DIR}/SNPsplit_reports
    # Parental genomes
    mkdir -p ${OUT_DIR}/parental_genomes
    mv *_${STRAIN}_*.gz ${STRAIN}.fa ${OUT_DIR}/parental_genomes
    # Delete files
    rm -r SNPs_${STRAIN}/ ${STRAIN}_full_sequence
done
