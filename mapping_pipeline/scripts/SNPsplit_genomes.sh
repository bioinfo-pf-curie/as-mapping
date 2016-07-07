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

# -- Run SNPsplit_genome_preparation to generate parental and N-masked genomes
echo "Generating parental and N-masked genomes ..."

${SNPsplit_gen} --nmasking --strain ${id_geno1} --strain2 ${id_geno2} --reference_genome ${ref_dir} --vcf_file ${full_vcf}

# -- Concatenetion of the parental genome in fasta files
for geno in $id_geno1 $id_geno2
do
	for i in `seq 1 19` X Y M 
	do
		sed 's/>/>chr/' ${geno}_full_sequence/chr${i}.SNPs_introduced.fa >> ${fasta_out}${geno}.fa
	done
done

# -- Concatenation of the N-masked genome
for i in `seq 1 19` X Y M
do
	sed 's/>/>chr/' ${id_geno1}_${id_geno2}_*_N-masked/chr${i}.N-masked.fa >> ${fasta_out}N-masked_${id_geno1}_${id_geno2}.fa
done

# -- Deleting non used files
mkdir reports
mv *.txt reports/
rm -r *_full_sequence/ *_N-masked/ SNPs_*/ *.gz

# -- Generation of bowtie2 indexes
echo "Creating bowtie2 indexes ..."

${bowtie2}bowtie2-build -f ${fasta_out}${id_geno1}.fa ${bowtie2_indexes}${id_geno1}
${bowtie2}bowtie2-build -f ${fasta_out}${id_geno2}.fa ${bowtie2_indexes}${id_geno2}
${bowtie2}bowtie2-build -f ${fasta_out}N-masked_${id_geno1}_${id_geno2}.fa ${bowtie2_indexes}N-masked_${id_geno1}_${id_geno2}


end=`date +%s`
echo " ==================================== "
echo "Script run in "$((end-start))" seconds."
echo " ==================================== "
