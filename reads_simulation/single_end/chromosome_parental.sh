#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#	chromosome_parental.sh : Script to generate the fasta sequence of a chromosome for a given strain (here chrX)

#### Parameters #### --------------------------------------------------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-i) id_geno=$2; shift;;
    (-v) vcf_geno=$2; shift;;
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

# Temporary file
tmp=${RANDOM}.tmp

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-i"" <name of the strain>"
    echo -e "-v"" <VCF/BED file with SNPs of the strain>"
    echo -e "-c"" <Config file>"
	echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------------------------------------------------

mkdir -p $fasta_outdir $vcf_outdir

if [[ ! -e ${vcf_geno} ]]; then echo "Error: VCF file not found. Exit"; exit; fi
if [[ -z ${id_geno} ]]; then echo "Error : Please specify the name of the strain. Exit"; exit; fi

# Select only SNPs that are homozygous in the vcf file
homo_vcf=${vcf_outdir}${chr}_homo_$(basename ${vcf_geno})
if [[ ! -e ${homo_vcf} ]]
then
	echo -e "  |\t$(basename $0) : Selecting homozygous SNPs of ${id_geno} from VCF on chr${chr}..."
	awk -v chr=${chr} '{if ($1 ~ /^#/) print; else if ($1 == chr) {split ($10,a,":"); if (a[1]=="1/1") print}}' ${vcf_geno} > ${homo_vcf}
fi

# Run vcf2diploid to generate parental genomes
echo -e "  |\t$(basename $0) : Generating parental chromosome of ${id_geno} on chr${chr}..."
java -jar $vcf2diploid -id ${id_geno} -chr ${ref} -vcf ${homo_vcf} -outDir $fasta_outdir

# Cleaning (Only paternal.fa is kept as both generated genomes should be identical due to the presence of only homozygous SNPs)
rm ${fasta_outdir}chr${chr}_${id_geno}_maternal.fa ${fasta_outdir}chr${chr}_${id_geno}.map ${fasta_outdir}*.chain

# Rename paternal genome as the alternatif genome
mv ${fasta_outdir}chr${chr}_${id_geno}_paternal.fa ${fasta_outdir}chr${chr}_${id_geno}.fa
sed "s/_paternal//" ${fasta_outdir}chr${chr}_${id_geno}.fa > ${tmp} && mv ${tmp} ${fasta_outdir}chr${chr}_${id_geno}.fa
