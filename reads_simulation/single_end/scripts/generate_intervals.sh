#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#   generate_intervals.sh : This script creates intervals around SNPs of a chosen size in BED format :
#       - From a VCF file ($diff_vcf)
#       - Select only intervals with mappability of 1
#       - Only in non intergenic regions

#### Parameters #### ---------------------------------------------------------------------------------------------------------------

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

# Size of intervals used for the simulation
inter_size=$read_length # Here the interval chosen depends on the read length chosen for the ART simulation

# Temporary directory and files 
tmp_dir=${tmp_outdir}/$RANDOM
tmp_bed=${tmp_dir}/intervals.bed.tmp
map_uniq=${tmp_dir}/uniq.map.tmp

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}


#### Main #### --------------------------------------------------------------------------------------------------------------------

mkdir -p ${vcf_outdir} ${tmp_outdir} ${tmp_dir}

# Generation of VCF file with different SNPs between the two strain if it does not already exists
vcf_name=$(basename $full_vcf)
diff_vcf=${vcf_outdir}/${vcf_name%.vcf}_${id_geno1}_${id_geno2}.vcf
if [[ ! -e $diff_vcf ]]
then
	echo -e "  |\t$(basename $0) : Generating VCF file of different SNPs between $id_geno1 and $id_geno2 ..."
	$extract_SNPs -i $full_vcf -r $id_geno1 -a $id_geno2 -f 1 > $diff_vcf
fi

# BED file generation
echo -e "  |\t$(basename $0) : Generating intervals BED file ..."

# Awk to create intervals of 2*$read_length around every SNPs of the VCF on the desired chromosome
awk -v i=$inter_size '{if ($1 !~ /^#/) {print "chr"$1,$2-i,$2+(i-1),$3":"$4"/"$5}}' OFS='\t' ${diff_vcf} > ${tmp_bed}

# Selection of intervals on mappability if bed specified by user
if [[ -e ${mappa} ]]
then
	# Selection of region with mappability == 1
	awk '{if ($4==1) print}' OFS='\t' $mappa > $map_uniq
	${bedtools} intersect -a ${tmp_bed} -b ${map_uniq} -wa -f 1 > ${tmp_dir}/tmp
	mv ${tmp_dir}/tmp ${tmp_bed}
	# -f 1 means the mappability has to be 1 on all the interval
fi

# Selection of regions of interest (specified in bed file)
${bedtools} intersect -a ${tmp_bed} -b ${regions} -wa -f 1 > ${tmp_dir}/tmp
mv ${tmp_dir}/tmp ${tmp_bed}

# Split interval BED files for allele specific read simulation
echo -e "  |\t$(basename $0) : Defining ratio for intervals ..."

if [[ -e ${ASratio} ]]
then
	while read line
	do
		ratio=`echo "$line" | cut -f5`
    	echo "$line" >> ${tmp_dir}/${ratio}.ratio
	done < ${ASratio}
	for r in `ls --color=never ${tmp_dir}/*.ratio`
	do
		$bedtools intersect -a ${tmp_bed} -b ${r} -wa -f 1 > ${r%.ratio}.bed
	done
	$bedtools intersect -a ${tmp_bed} -b ${tmp_dir}/*.ratio -wa -v >> ${tmp_dir}/0.5.bed
else
	mv ${tmp_bed} ${tmp_dir}/0.5.bed
fi

# Generate reads
echo -e "  |\t$(basename $0) : Generating reads with ART ..."

for bed in `ls --color=never ${tmp_dir}/*.bed`
do
	nb_geno1=$(echo "${coverage} * $(basename ${bed%.bed})" | bc)
	nb_geno2=$(echo "${coverage} * (1-$(basename ${bed%.bed}))" | bc)
	${simreads}/scripts/generate_reads.sh -i ${id_geno1} -b ${bed} -n ${nb_geno1} -c ${config}
	${simreads}/scripts/generate_reads.sh -i ${id_geno2} -b ${bed} -n ${nb_geno2} -c ${config}
done

# Cleaning the file
rm -r ${tmp_dir}
