#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :

#### Parameters #### --------------------------------------------------------------------------------------------------------------

while [ $# -gt 0 ] 
do
    case "$1" in
    (-c) config=$2; shift;;
	(-n) number_reads=$2; shift;;
    (-h) usage;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;; 
    (*)  break;;
    esac
    shift
done

if [[ -z ${config} || ! -e ${config} ]]; then echo "Error : Please specify a config file. Exit"; exit; fi
source ${config}

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
	echo -e "-n"" <Number of reads to be generated from each allele>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------------------------------------------------

start=`date +%s`

mkdir -p ${main_out}

#### STEP 1 : Generation of the parental chromosomes

if [[ ! -e ${fasta_outdir}chr${chr}_${id_geno1}.fa ]]
then
	echo "Generation of the parental chromosome for ${id_geno1} ..."
	${simreads}scripts/build_parental_reference.sh -i ${id_geno1} -c ${config}
fi
if [[ ! -e ${fasta_outdir}chr${chr}_${id_geno2}.fa ]]
then
	echo "Generation of the parental chromosome for ${id_geno2} ..."
	${simreads}scripts/build_parental_reference.sh -i ${id_geno2} -c ${config}
fi


#### STEP 2 : Generate intervals in BED format (non intergenic region and mappabality == 1 on reference genome)
if [[ ! -e ${bed_outdir}${chr}_${int_bed} ]]
then
	echo "Generation of the intervals in BED format ..."
	mkdir -p ${bed_outdir}
	${simreads}scripts/generate_intervals.sh -c ${config}
fi


#### STEP 3 : Generate a chosen number of reads on each interval for the two strains
echo "Generation of reads with ART ..."
${simreads}scripts/generate_reads.sh -i ${id_geno1} -b ${bed_outdir}${chr}_${int_bed} -n ${number_reads} -c ${config}
${simreads}scripts/generate_reads.sh -i ${id_geno2} -b ${bed_outdir}${chr}_${int_bed} -n ${number_reads} -c ${config}


#### STEP 4 : Cleaning of the different outputs (rephasing, removing unused files, merging .fq and .sam)
echo "Merging and Cleaning files ..."
${simreads}scripts/cleaning_outputs.sh -c ${config}

end=`date +%s`
echo " ==================================== "
echo "Script run in "$((end-start))" seconds."
echo " ==================================== "
