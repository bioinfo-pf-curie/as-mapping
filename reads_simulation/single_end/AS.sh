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
	echo -e "-n"" <Number of reads to be generated>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------------------------------------------------

start=`date +%s`

mkdir -p ${main_out}

#### STEP 1 : Generation of the parental X chromosomes
if [[ ! -e ${fasta_outdir}${id_geno1}.fa ]]
then
	echo "Generation of the parental chromosome for ${id_geno1} ..."
	${simreads}build_parental_reference.sh -v ${vcf_geno1} -i ${id_geno1} -c ${config}
fi
if [[ ! -e ${fasta_outdir}${id_geno2}.fa ]]
then
	echo "Generation of the parental chromosome for ${id_geno2} ..."
	${simreads}build_parental_reference.sh -v ${vcf_geno2} -i ${id_geno2} -c ${config}
fi


#### STEP 2 : Generate intervals in BED format (non intergenic region and mappabality == 1 on reference genome)
echo "Generation of the intervals in BED format ..."
mkdir -p ${bed_outdir}
if [[ ! -e ${bed_outdir}${int_bed} ]]
then
	${simreads}generate_intervals.sh -c ${config}
fi

# optional STEP (to create AS) : Separate interval BED file into several BED files

if [[ ! -e ${bed_outdir}other_${int_bed} || ! -e ${bed_outdir}${id_geno1}${int_bed} || ! -e ${bed_outdir}${id_geno2}${int_bed} ]]
then
	${bedtools} intersect -a ${bed_outdir}${chr}_${int_bed} -b ${ASspecgeno1} -wa | uniq > ${bed_outdir}${id_geno1}spec_${chr}_${int_bed}
	${bedtools} intersect -a ${bed_outdir}${chr}_${int_bed} -b ${ASspecgeno2} -wa | uniq > ${bed_outdir}${id_geno2}spec_${chr}_${int_bed}
	${bedtools} intersect -a ${bed_outdir}${chr}_${int_bed} -b ${ASspecgeno1} ${ASspecgeno2} -v | uniq > ${bed_outdir}other_${chr}_${int_bed}
fi

#### STEP 3 : Generate a chosen number of reads on each interval for the two strains
echo "Generation of reads with ART ..."
# AS geno1
${simreads}generate_reads.sh -i ${id_geno1} -b ${bed_outdir}${id_geno1}spec_${chr}_${int_bed} -n ${number_reads} -c ${config}
${simreads}generate_reads.sh -i ${id_geno2} -b ${bed_outdir}${id_geno1}spec_${chr}_${int_bed} -n $((number_reads/10)) -c ${config}
# AS geno2
${simreads}generate_reads.sh -i ${id_geno1} -b ${bed_outdir}${id_geno2}spec_${chr}_${int_bed} -n $((number_reads/10)) -c ${config}
${simreads}generate_reads.sh -i ${id_geno2} -b ${bed_outdir}${id_geno2}spec_${chr}_${int_bed} -n ${number_reads} -c ${config}
# no AS
${simreads}generate_reads.sh -i ${id_geno1} -b ${bed_outdir}other_${chr}_${int_bed} -n $((number_reads/2)) -c ${config}
${simreads}generate_reads.sh -i ${id_geno2} -b ${bed_outdir}other_${chr}_${int_bed} -n $((number_reads/2)) -c ${config}

#### STEP 4 : Cleaning of the different outputs (rephasing, removing unused files, merging .fq and .sam)
echo "Merging and Cleaning files ..."
${simreads}cleaning_outputs.sh -c ${config}

end=`date +%s`
echo " ==================================== "
echo "Script run in "$((end-start))" seconds."
echo " ==================================== "
