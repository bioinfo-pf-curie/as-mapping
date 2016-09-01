#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to map reads on reference genome using Bowtie2
#		STEP 1 : Mapping of the reads on reference genome
#		STEP 2 : BAM processing for Allele Specific analysis


#### Parameters #### --------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-c) config=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;; 
    (*)  break;;
    esac
    shift
done

if [[ -z $config ]]
then
	echo "ERROR : you need to specify a config file. Exit."
	exit
fi
source ${config}

bam_analysis=${map_path}scripts/ref/ref_analysis.sh

#### Function #### ----------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------

start=`date +%s`

# Set up output directory for this method of mapping in the main output directory
BAM_OUT=${OUT_DIR}/mapping_reference
mkdir -p ${BAM_OUT}/

# Get the name of the reference genome to use already existing indexes (might change this part to work directly on my own indexes)
ID_REF=$(basename ${REF_GENO})
ID_REF=${ID_REF%.fa}

# Create bowtie2 indexes if they do not exist
if [ ! -e ${INDEXES}${ID_REF}.rev.2.bt2 ]
then
	echo "Generating bowtie2 indexes ..."
	mkdir $INDEXES
	${bowtie2}bowtie2-build -f ${ref_geno} ${INDEXES}${ID_REF}
fi

${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -x ${INDEXES}${ID_REF} -U $fq_reads | ${samtools} view -bS - > ${sam_out}${ID_REF}.bam


##### STEP 2 : BAM analysis -------------------------------------------------------------

${bam_analysis} -c ${config}

end=`date +%s`
echo "## mapping_reference.sh run in "$((end-start))" seconds."
