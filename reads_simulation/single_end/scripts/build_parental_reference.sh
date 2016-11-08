#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#	build_parental_reference.sh : Script to generate the fasta sequence of a chromosome for a given strain

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-i"" <name of the strain>"
    echo -e "-c"" <Config file>"
	echo -e "-h"" <help>"
    exit
}

#### Parameters #### --------------------------------------------------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-i) id_geno=$2; shift;;
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


#### Main #### --------------------------------------------------------------------------------------------------------------------

mkdir -p ${main_out} ${fasta_outdir}
 
if [[ -z ${id_geno} ]]; then echo "Error : Please specify the name of the strain. Exit"; exit; fi


if [[ ${id_geno} == 'REF' ]]
then
    for i in `ls -d --color=never ${ref_dir}/*`
    do  
        if [[ $i =~ MT ]]
        then
            sed 's/>MT/>chrM/' ${i} >> ${fasta_outdir}/${id_geno}.fa
        else
            sed 's/>/>chr/' ${i} >> ${fasta_outdir}/${id_geno}.fa
        fi
    done
else
    # Run SNPsplit to generate parental genome
    ${SNPsplit_gen} --no_nmasking --strain ${id_geno} --reference_genome ${ref_dir} --vcf_file ${full_vcf} 

    for i in `ls -d --color=never ${id_geno}_full_sequence/*`
    do  
        if [[ $i =~ MT ]]
        then
            sed 's/>MT/>chrM/' ${i} >> ${fasta_outdir}/${id_geno}.fa
        else
            sed 's/>/>chr/' ${i} >> ${fasta_outdir}/${id_geno}.fa
        fi
    done
    # Cleaning (Only paternal.fa is kept as both generated genomes should be identical due to the presence of only homozygous SNPs)
    mkdir -p reports/
    mv *.txt reports/
    rm -r ${id_geno}_full_sequence/ SNPs_*/ *.gz
fi


