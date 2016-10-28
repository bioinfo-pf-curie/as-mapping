#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Script to build bowtie2 indexes


#### Function #### -------------------------------------------------------------------

function usage {
    echo -e "Usage : $0"
    echo -e "-f"" <FASTA file>"
    echo -e "-b"" <BOWTIE2 directory>"
    echo -e "-o"" <OUTPUT directory for indexes>"
    echo -e "-h"" <help>"
    exit
}

#### Parameters #### -----------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-f) FASTA=$2; shift;;
    (-b) BOWTIE2_DIR=$2; shift;;
    (-o) OUTDIR=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;; 
    (*)  break;;
    esac
    shift
done

if [[ -z $FASTA ]]
then
    echo "ERROR : you need to specify a FASTA file (-f). Exit."
    exit
fi

B2_BUILD_IND=${BOWTIE2_DIR}/bowtie2-build

#### Main #### -----------------------------------------------------------------------

mkdir -p $OUTDIR

NAME=$(basename ${FASTA%.*})
echo "$0: Building Bowtie2 indexes for ${NAME} in ${OUTDIR}"
${B2_BUILD_IND} -f ${FASTA} ${OUTDIR}/${NAME}
