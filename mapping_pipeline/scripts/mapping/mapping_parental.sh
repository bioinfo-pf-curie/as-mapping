#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to perform alignment to parental genomes
#		STEP 1 : Generation of parental genomes
#		STEP 2 : Mapping of the reads on both genomes
#		STEP 3 : Select best alignement between both genomes and BAM processing for Allele Specific analysis

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

bam_analysis=${map_path}scripts/par/par_analysis.sh

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
sam_out=$sam_out"mapping_parental/"
mkdir -p $sam_out $main_out


##### STEP 1 : Generation of parental genomes -------------------------------------------
#	Following this script, both parental genomes are generated in the $fasta_out directory
#	In parallel, N-masked genome is also generated

if [[ ! -e ${fasta_out}${id_geno1}.fa && ! -e ${fasta_out}${id_geno2}.fa ]]
then # Create parental genomes and indexes (generation of indexes within the script)
	${SNPsplit_genomes} -c ${config}
fi


##### STEP 2 : Alignment with Bowtie2 ---------------------------------------------------

${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -x $bowtie2_indexes$id_geno1 -U $fq_reads | ${samtools} view -bS - > ${sam_out}${id_geno1}.bam
${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -x $bowtie2_indexes$id_geno2 -U $fq_reads | ${samtools} view -bS - > ${sam_out}${id_geno2}.bam


##### STEP 3 : Comparison to find best alignment and BAM analysis  ---------------------------

${bam_analysis} -c ${config}

end=`date +%s`
echo " ====================================================== "
echo "  mapping_parental.sh run in "$((end-start))" seconds."
echo " ====================================================== "
