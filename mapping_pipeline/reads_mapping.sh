#!/usr/bin/env bash
# Author(s) : Kenzo-Hugo Hillion
# Comment(s) : Main part of the mapping pipeline

#### Function #### ----------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <CONFIG file>"
    echo -e "-f"" <Forward reads (paired-end) or reads (single-end)>"
    echo -e "-r"" <[Paired-end ONLY] reverse reads>"
    echo -e "-o"" <Output directory for the alignment files>"
    echo -e "-n"" <Output name for the alignment files>"
    echo -e "-h"" <help>"
    exit
}

#### Parameters #### --------------------------------------------------------------------

while [[ $# -gt 0 ]]
do
    case "$1" in
    (-c) CONFIG=$2; shift;;
    (-f) FQ_READS_F=$2; shift;;
    (-r) FQ_READS_R=$2; shift;;
    (-o) OUT_DIR=$2; shift;;
    (-n) OUT_NAME=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: ERROR - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

if [[ -z $CONFIG ]]; then echo "$0: ERROR - you need to specify a CONFIG file (-c). Exit." 1>&2; exit 1; fi
source $CONFIG
source ${PIPELINE_PATH}/includes/path_fct.inc

#### Main #### --------------------------------------------------------------------------

# Create output directory
mkdir -p ${OUT_DIR}/
# Run mapping strategies

if [[ ${MAP_REF} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_REFERENCE_B2} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_REFERENCE_B2} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_REFERENCE_TOPHAT} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_REFERENCE_TOPHAT} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to reference genome SKIPPED." 1>&2
    fi
fi

if [[ ${MAP_N} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_NMASKED_B2} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
	    else
            ${MAPPING_NMASKED_B2} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_NMASKED_TOPHAT} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_NMASKED_TOPHAT} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to N-masked genome SKIPPED." 1>&2
    fi
fi

if [[ ${MAP_PAR} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_PARENTAL_B2} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_PARENTAL_B2} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_PARENTAL_TOPHAT} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_PARENTAL_TOPHAT} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to parental genomes SKIPPED." 1>&2
    fi
fi

if [[ ${MAP_DIP} -eq 1 ]]
then
    if [[ ${MAPPER} == 'BOWTIE2' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_DIPLOID_B2} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_DIPLOID_B2} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    elif [[ ${MAPPER} == 'TOPHAT' ]]
    then
        if [[ -e ${FQ_READS_R} ]]
        then
            ${MAPPING_DIPLOID_TOPHAT} -f $FQ_READS_F -r $FQ_READS_R -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        else
            ${MAPPING_DIPLOID_TOPHAT} -f $FQ_READS_F -o $OUT_DIR -n $OUT_NAME -c $CONFIG
        fi
    else
        echo "$0: WARNING - The selected mapper [${MAPPER}] does not exist." 1>&2
        echo "$0: WARNING - Mapping to diploid genome SKIPPED." 1>&2
    fi
fi

exit 0