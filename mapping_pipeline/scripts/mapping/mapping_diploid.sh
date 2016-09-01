#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to perform alignment to a diploid genome
#		STEP 1 : Generation of the diploid genome
#		STEP 2 : Mapping of the reads on the genome
#		STEP 3 : BAM processing for Allele Specific analysis

#### Parameters #### --------------------------------------------------------------------

# -- Local scripts for BAM analysis
selectBest=${map_path}scripts/par/select_best.py

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

bam_analysis=${map_path}scripts/dip/dip_analysis.sh

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
sam_out=$sam_out"mapping_diploid/"
mkdir -p $sam_out


##### STEP 1 : Generation of the diploid genome -----------------------------------------
#	Following this script, both parental genomes are generated in the $fasta_out directory
#	with the named specified in $fasta_geno1 and $fasta_geno2

if [[ ! -e ${fasta_out}${id_geno1}_${id_geno2}.fa ]]
then # Create diploid genome and indexes (generation of indexes within the script)
	${diploid_genome} -c ${config}
fi


##### STEP 2 : Alignment with Bowtie2 --------------------------------------------------------

${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -k 3 -x $bowtie2_indexes${id_geno1}_${id_geno2} -U $fq_reads | ${samtools} view -bS - > ${sam_out}${id_geno1}_${id_geno2}.bam


##### STEP 3 : BAM analysis TO BE DONE  -------------------------------------------------

${bam_analysis} -c ${config}

end=`date +%s`
echo " ====================================================== "
echo "  mapping_diploid.sh run in "$((end-start))" seconds."
echo " ====================================================== "
