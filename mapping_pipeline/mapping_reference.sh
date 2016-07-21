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

##### STEP 1 : Alignment to reference genome --------------------------------------------

# Set up output directory for this method of mapping in the main output directory
sam_out=${sam_out}mapping_reference/
mkdir -p ${sam_out}

# Get the name of the reference genome to use already existing indexes (might change this part to work directly on my own indexes)
id_ref=$(basename ${ref_geno})
id_ref=${id_ref%.fa}

# Create bowtie2 indexes if they do not exist
if [ ! -e ${bowtie2_indexes}${id_ref}.rev.2.bt2 ]
then
	echo "Generating bowtie2 indexes ..."
	mkdir $bowtie2_indexes
	${bowtie2}bowtie2-build -f ${ref_geno} ${bowtie2_indexes}${id_ref}
fi

${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -x ${bowtie2_indexes}${id_ref} -U $fq_reads | ${samtools} view -bS - > ${sam_out}${id_ref}.bam


##### STEP 2 : BAM analysis -------------------------------------------------------------

${bam_analysis} -c ${config}

end=`date +%s`
echo "## mapping_reference.sh run in "$((end-start))" seconds."
