#!/usr/bin/env bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Script to reconstruct parental genome in fasta format from VCF files containing SNPs

#### Functions #### ------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-s/--strain"" <Name of the strain>" 
    echo -e "-v/--vcf"" <VCF file>" 
    echo -e "-g/--snpsplit_gen"" <SNPsplit_genome_preparation path>"
    echo -e "-r/--reference"" <Directory with reference chromosomes (Ensembl format)>" 
    echo -e "-o/--outdir"" <Output directory>" 
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
        --strain) args="${args}-s ";;
        --vcf) args="${args}-v ";;
        --snpsplit_gen) args="${args}-g ";;
        --reference) args="${args}-r ";;
        --outdir) args="${args}-o ";;
        *) [[ "${arg:0:1}" == "-" ]] || delim="\""
            args="${args}${delim}${arg}${delim} ";;
    esac
done

eval set -- ${args}

# Short arguments parsing
while getopts :s:v:g:r:o:h option
do
    case "${option}" in
        s) STRAIN=${OPTARG};;
        v) VCF=${OPTARG};;
        g) SNPSPLIT_GEN=${OPTARG};;
        r) REF_DIR=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        h) usage;;
        \?) echo "Invalid option: - ${OPTARG}"; exit 1;;
        :) echo "Option -${OPTARG} requires an argument." >&2; exit 1;;
    esac
done


#### Main #### -----------------------------------------------------------------------

mkdir -p $OUT_DIR

# -- Run SNPsplit_genome_preparation to generate parental and N-masked genomes
echo "$0: Generating parental genome of ${STRAIN} ..."
echo "in $OUT_DIR"

${SNPSPLIT_GEN} --strain ${STRAIN} --reference_genome ${REF_DIR} --vcf_file ${VCF} --no_nmasking

# -- Concatenetion of the parental genome in fasta files
for i in `ls -d --color=never ${STRAIN}_full_sequence/*`
do
    if [[ $i =~ MT ]]
    then
        sed 's/>MT/>chrM/' ${i} >> ${OUT_DIR}/${STRAIN}.fa
    else
        sed 's/>/>chr/' ${i} >> ${OUT_DIR}/${STRAIN}.fa
	fi
done

# -- Save and delete files
# Reports
mkdir -p ${OUT_DIR}/SNPsplit_reports
mv ${STRAIN}_*report.txt ${OUT_DIR}/SNPsplit_reports
# SNPs
mv *_${STRAIN}_*.gz ${OUT_DIR}
# Delete files
rm ${STRAIN}_specific_*.txt
rm -r SNPs_${STRAIN} ${STRAIN}_full_sequence

exit 0
