#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :

#### Parameters #### --------------------------------------------------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-i) id_geno=$2; shift;;
    (-b) bed=$2; shift;;
	(-n) reads_number=$2; shift;;
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

# temporary files
tmp_int=$RANDOM.inttmp

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-b"" <BED file of intervals>"
    echo -e "-i"" <Name of the strain>"
	echo -e "-n"" <Number of reads to be generated for each intervals>"
    echo -e "-c"" <Config file>"
	echo -e "-h"" <help>"
    exit
}


#### Main #### --------------------------------------------------------------------------------------------------------------------


mkdir -p ${art_outdir} ${bed_outdir}

# Generation of the different fragments in fasta around the SNPs
fasta=chr${chr}_$id_geno".fa"
int_fa=${id_geno}_$(basename $bed)
int_fa=${int_fa%.bed}.fa

if [[ ! -e ${fasta_outdir}${int_fa} ]]
then
	echo -e "  |\t$(basename $0) : Getting fasta from interval bed file ..."
	${bedtools} getfasta -fi ${fasta_outdir}${fasta} -bed ${bed} -fo ${fasta_outdir}${int_fa}
	awk -v id=${id_geno} '{if ($1 ~ /^>/) print $1"_"id; else print}' ${fasta_outdir}${int_fa} > ${fasta_outdir}${tmp_int} && mv ${fasta_outdir}${tmp_int} ${fasta_outdir}${int_fa}
fi

# Simulation of reads with ART
echo -e "  |\t$(basename $0) : Generation of reads with ART ..."
${art} -i ${fasta_outdir}${int_fa} -o ${art_outdir}${int_fa%.fa} -l ${read_length} -c ${reads_number} -ss ${sequencer} -rs ${rs} -ir ${ir} -dr ${dr} -qs ${subs} -sam
