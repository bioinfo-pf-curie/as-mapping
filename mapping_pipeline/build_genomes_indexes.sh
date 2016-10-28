#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Build necessary files from a config file

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-h"" <help>"
    exit
}


#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-c) CONFIG=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

if [[ -z $CONFIG ]]; then echo "ERROR : you need to specify a CONFIG file (-c). Exit."; exit 1; fi
source $CONFIG
source ${PIPELINE_PATH}/includes/path_fct.inc


#### Main #### --------------------------------------------------------------------------

# Build missing files for each strategies

if [[ ${MAP_REF} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]] || [[ ${MAPPER} == 'TOPHAT' ]] # Bowtie2 or Tophat
    then
        if [[ -z ${B2_INDEX_REF} ]]
        then # Make Bowtie2 indexes if user did not specify path
            ${BUILD_B2_INDEX} -f ${REF_GENO} -b ${BOWTIE2_DIR} -o ${INDEXES}
        fi
    fi
fi

if [[ ${MAP_N} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]] || [[ ${MAPPER} == 'TOPHAT' ]] # Bowtie2 or Tophat
    then
        if [[ -z ${B2_INDEX_NMASK} ]]
        then # Make Bowtie2 indexes if user did not specify path
            if [[ -z ${FASTA_NMASK} ]]
            then # Build N-masked genome first
                ${BUILD_NMASK_GENOME} -p ${ID_GENO1} -m ${ID_GENO2} -v ${FULL_VCF} -g ${SNPSPLIT_GEN} -r ${REF_DIR} -o ${FASTA_OUT}
                FASTA_NMASK=${FASTA_OUT}/N-masked_${ID_GENO1}_${ID_GENO2}.fa
            fi
            ${BUILD_B2_INDEX} -f ${FASTA_NMASK} -b ${BOWTIE2_DIR} -o ${INDEXES} 
        fi
    fi
fi

if [[ ${MAP_DIP} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]] || [[ ${MAPPER} == 'TOPHAT' ]] # Bowtie2 or Tophat
    then
        # Diploid
        if [[ -z ${B2_INDEX_DIP} ]]
        then
            if [[ -z ${FASTA_DIP} ]]
            then
                ${BUILD_DIPLOID_GENOME} -p ${ID_GENO1} -m ${ID_GENO2} -v ${FULL_VCF} -g ${SNPSPLIT_GEN} -r ${REF_DIR} -o ${FASTA_OUT}
                FASTA_DIP=${FASTA_OUT}/${ID_GENO1}_${ID_GENO2}.fa
                FASTA_GENO1=${FASTA_OUT}/${ID_GENO1}.fa
                FASTA_GENO2=${FASTA_OUT}/${ID_GENO2}.fa
            fi
            ${BUILD_B2_INDEX} -f ${FASTA_DIP} -b ${BOWTIE2_DIR} -o ${INDEXES}
        fi 
    fi
fi


if [[ ${ID_GENO1} == 'REF' ]]
then
    B2_INDEX_GENO1=${B2_INDEX_REF}
    FASTA_GENO1=${REF_GENO}
elif [[ ${ID_GENO2} == 'REF' ]]
then
    B2_INDEX_GENO2=${B2_INDEX_REF}
    FASTA_GENO2=${REF_GENO}
fi
if [[ ${MAP_PAR} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]] || [[ ${MAPPER} == 'TOPHAT' ]] # Bowtie2 or Tophat
    then
        # First genotype
        if [[ -z ${B2_INDEX_GENO1} ]]
        then # Make Bowtie2 indexes for GENO1
            if [[ -z ${FASTA_GENO1} ]]
            then
                ${BUILD_PARENTAL_GENOME} -s ${ID_GENO1} -v ${FULL_VCF} -g ${SNPSPLIT_GEN} -r ${REF_DIR} -o ${FASTA_OUT}
                FASTA_GENO1=${FASTA_OUT}/${ID_GENO1}.fa
            fi
            ${BUILD_B2_INDEX} -f ${FASTA_GENO1} -b ${BOWTIE2_DIR} -o ${INDEXES}
        fi 
        # Second genotype
        if [[ -z ${B2_INDEX_GENO2} ]]
        then # Make Bowtie2 indexes for GENO2
            if [[ -z ${FASTA_GENO2} ]]
            then
                ${BUILD_PARENTAL_GENOME} -s ${ID_GENO2} -v ${FULL_VCF} -g ${SNPSPLIT_GEN} -r ${REF_DIR} -o ${FASTA_OUT}
                FASTA_GENO2=${FASTA_OUT}/${ID_GENO2}.fa
            fi
            ${BUILD_B2_INDEX} -f ${FASTA_GENO2} -b ${BOWTIE2_DIR} -o ${INDEXES}
        fi
    fi
fi

exit 0
