#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : Script to reconstruct two parental genomes in fasta format from VCF files containing SNPs

#### Parameters #### --------------------------------------------------------------------------------------------------------------

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

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------------------------------------------------

start=`date +%s`
mkdir -p $fasta_out $fasta_out$bowtie2_indexes

# Generation of parental genomes
if [[ ! -e $fasta_out$fasta_geno1 && ! -e $fasta_out$fasta_geno2 ]]
then # Create parental genomes and indexes (generation of indexes within the script)
    ${parental_genomes} -c ${config}
fi

# Renaming chromosomes to correspond to its origin and concatenate to form the diploid genome
fasta_diploid=${id_geno1}_${id_geno2}.fa && rm ${fasta_out}${fasta_diploid}
for id in ${id_geno1} ${id_geno2}
do
	awk -v id=$id  '{if ($1 ~ /^>/) print $1"_"id; else print}' $fasta_out${id}.fa >> ${fasta_out}${fasta_diploid}
done

# -- Generation of bowtie2 indexes
echo "Creating bowtie2 indexes ..."
bowtie2_indexes=$fasta_out$bowtie2_indexes
${bowtie2}bowtie2-build -f ${fasta_out}${fasta_diploid} ${bowtie2_indexes}${id_geno1}_${id_geno2}


end=`date +%s`
echo " =================================================== "
echo "  diploid_genome.sh run in "$((end-start))" seconds."
echo " =================================================== "
