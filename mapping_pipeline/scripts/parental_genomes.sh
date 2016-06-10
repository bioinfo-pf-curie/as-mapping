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
mkdir -p $fasta_out $vcf_out $fasta_out$bowtie2_indexes


# -- Extraction of only homozygous SNPs with awk
echo "Selecting homozygous SNPs from VCF ..."

if [[ ! -e $vcf_out$homo_vcf1 ]]
then # Create the file only if it is not present already
	awk '{if ($1 ~ /^#/) print; else{split ($10,a,":"); if (a[1]=="1/1") print}}' $vcf_geno1 > $vcf_out$homo_vcf1
fi
if [[ ! -e $vcf_out$homo_vcf2 ]] 
then # Create the file only if it is not present already
	awk '{if ($1 ~ /^#/) print; else{split ($10,a,":"); if (a[1]=="1/1") print}}' $vcf_geno2 > $vcf_out$homo_vcf2
fi

# -- vcf2diploid applied on the reference genome using vcf files (with homozygous SNPs of the strain)
echo "Generating parental genomes ..."

mkdir ${fasta_out}chromosomes$id_geno1
mkdir ${fasta_out}chromosomes$id_geno2
java -Xmx40G -jar $vcf2diploid -id $id_geno1 -chr $ref_geno -vcf $vcf_out$homo_vcf1 -outDir ${fasta_out}chromosomes$id_geno1
java -Xmx40G -jar $vcf2diploid -id $id_geno2 -chr $ref_geno -vcf $vcf_out$homo_vcf2 -outDir ${fasta_out}chromosomes$id_geno2

# -- Concatenetion of the genome in one fasta file (from paternal one, which is supposed to be identical from maternal in this case)
for fasta in $fasta_geno1 $fasta_geno2
do
	for i in `seq 1 19` X Y M 
	do
		id=${fasta%.fa}
		sed -i 's/\_paternal//' ${fasta_out}chromosomes$id/chr${i}_${id}_paternal.fa 
		cat ${fasta_out}chromosomes$id/chr${i}_${id}_paternal.fa >> ${fasta_out}${fasta}
		# Renaming each chromosome by removing "_paternal"
		mv ${fasta_out}chromosomes$id/chr${i}_${id}_paternal.fa ${fasta_out}chromosomes$id/chr${i}_${id}.fa
	done
done

# -- Deleting non used files
for id in $id_geno1 $id_geno2
do
	rm ${fasta_out}chromosomes$id/*maternal.fa
	rm ${fasta_out}chromosomes$id/*.map
	rm ${fasta_out}chromosomes$id/*.chain
done

# -- Generation of bowtie2 indexes
echo "Creating bowtie2 indexes ..."
bowtie2_indexes=$fasta_out$bowtie2_indexes

${bowtie2}bowtie2-build -f $fasta_out$fasta_geno1 $bowtie2_indexes$id_geno1
${bowtie2}bowtie2-build -f $fasta_out$fasta_geno2 $bowtie2_indexes$id_geno2


end=`date +%s`
echo " ==================================== "
echo "Script run in "$((end-start))" seconds."
echo " ==================================== "
