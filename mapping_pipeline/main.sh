#!/usr/bin/env bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Main part of the mapping pipeline

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
    if [[ ${MAPPER} -eq 1 ]] # Bowtie2
    then
        ${MAPPING_REFERENCE_B2} -c $CONFIG
    elif [[ ${MAPPER} -eq 2 ]] # Tophat
    then
        ${MAPPING_REFERENCE_TOPHAT} -c $CONFIG
    fi
fi

if [[ ${MAP_N} -eq 1 ]]
then
    if [[ ${MAPPER} -eq 1 ]] # Bowtie2
    then
        ${MAPPING_NMASKED_B2} -c $CONFIG
    elif [[ ${MAPPER} -eq 2 ]] # Tophat
    then
        ${MAPPING_NMASKED_TOPHAT} -c $CONFIG
    fi
fi

if [[ ${MAP_PAR} -eq 1 ]]; then ${MAPPING_PARENTAL} -c $CONFIG ; fi
if [[ ${MAP_DIP} -eq 1 ]]; then ${MAPPING_DIPLOID} -c $CONFIG ; fi
