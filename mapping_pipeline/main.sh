#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Main part of the mapping pipeline

#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-c) CONFIG=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: ERROR - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

if [[ -z $CONFIG ]]; then echo "$0: ERROR - you need to specify a CONFIG file (-c). Exit." 1>&2; exit 1; fi
source $CONFIG

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------

# Create output directory
mkdir -p ${OUT_DIR}/
# Run mapping strategies

if [[ ${MAP_REF} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        ${MAPPING_REFERENCE_B2} -c $CONFIG
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        ${MAPPING_REFERENCE_TOPHAT} -c $CONFIG
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to reference genome SKIPPED." 1>&2
    fi
fi

if [[ ${MAP_N} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        ${MAPPING_NMASKED_B2} -c $CONFIG
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        ${MAPPING_NMASKED_TOPHAT} -c $CONFIG
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to N-masked genome SKIPPED." 1>&2
    fi
fi

if [[ ${MAP_PAR} -eq 1 ]]
then 
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        ${MAPPING_PARENTAL_B2} -c $CONFIG
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        ${MAPPING_PARENTAL_TOPHAT} -c $CONFIG
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to parental genomes SKIPPED." 1>&2
    fi
fi 

if [[ ${MAP_DIP} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        ${MAPPING_DIPLOID_B2} -c $CONFIG
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        ${MAPPING_DIPLOID_TOPHAT} -c $CONFIG
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to diploid genome SKIPPED." 1>&2
    fi
fi

exit 0 
