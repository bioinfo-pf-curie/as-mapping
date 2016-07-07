#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  		Script to perform alignment to a N masked genome
#		STEP 1 : Generation of the N-masked genome (masked bases are those that differ between the two studies strains
#		STEP 2 : Mapping of the reads on the N-masked genome
#		STEP 3 : BAM processing for Allele Specific studies

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

bam_analysis=${map_path}scripts/Nmask/N_analysis.sh

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
sam_out=$sam_out"mapping_N_masked/"
mkdir -p $sam_out

##### STEP 1 : Generation of N-masked genome --------------------------------------------

#       Following this script, N-masked genome is generated based on differential SNPs between 
#       the two strains specified by $id_geno1 and $id_geno2
#		In parallel, parental genomes are also generated

masked_genome="N-masked_"$id_geno1"_"$id_geno2
if [[ ! -e $fasta_out$masked_genome.fa ]]
then # Need to generate the masked genome
	${SNPsplit_genomes} -c ${config}
fi


##### STEP 2 : Alignment with Bowtie2 ---------------------------------------------------

${bowtie2}bowtie2 ${B2_OPTIONS} ${B2_SCORING_OPT} -x ${bowtie2_indexes}${masked_genome} -U $fq_reads | ${samtools} view -bS - > ${sam_out}${masked_genome}.bam


##### STEP 3 : BAM analysis -------------------------------------------------------------

${bam_analysis} -c ${config}

end=`date +%s`
echo " ======================================================= "
echo "  mapping_N_masked.sh run in "$((end-start))" seconds."
echo " ======================================================= "
