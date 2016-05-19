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

# Temporary files 
tmp_bed=$RANDOM.bedtmp
map_uniq=$RANDOM.maptmp

# Size of intervals used for the simulation
inter_size=$read_length # Here the interval chosen depends on the read length chosen for the ART simulation

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}


#### Main #### --------------------------------------------------------------------------------------------------------------------

# Generation of VCF file with different SNPs between the two strain if it does not already exists
if [[ -z $diff_vcf ]]
then
    vcf=`echo $full_vcf | tr '/' '\n' | tail -1`
    diff_vcf=$vcf_outdir${vcf%.vcf}$id_geno1"_"$id_geno2".vcf"
    if [[ ! -e $diff_vcf ]]
	then
    	echo -e "  |\t$(basename $0) : Generating VCF file of different SNPs between $id_geno1 and $id_geno2 ..."
		$extract_SNPs -i $full_vcf -r $id_geno1 -a $id_geno2 -f 1 > $diff_vcf
	fi
fi

# BED file generation
echo -e "  |\t$(basename $0) : Generating intervals BED file ..."
# Select non intergenic variants
vcf_non_inter=non_inter_$RANDOM.tmp
# Awk to create intervals of 2*100bp around every SNPs of a VCF files containing differential SNPs between 2 parental genomes in a BED file
awk -v chr=$chr '{if ($1 ~ /^#/) print; else if ($1 == chr) {split ($8,a,"|"); if (a[5] !~ /intergenic_variant/) print}}' OFS='\t' ${diff_vcf} > $vcf_non_inter

# Intervals from VCF file
awk -v i=$inter_size '{if ($1 !~ /^#/) print "chr"$1,$2-i,$2+(i-1),$3":"$4"/"$5}' OFS='\t' $vcf_non_inter > $tmp_bed 

# Selection of intervals on mappability if bed specified by user
if [[ -e ${mappa} ]]
then
	# Selection of region with mappability == 1
	awk '{if ($4==1) print}' OFS='\t' $mappa > $map_uniq
	${bedtools} intersect -a ${tmp_bed} -b ${map_uniq} -wa -f 1 > ${bed_outdir}${int_bed}
	# -f 1 means the mappability has to be 1 on all the interval
else
	mv $tmp_bed ${bed_outdir}${int_bed}
fi

# Cleaning the file
rm $tmp_bed $map_uniq $vcf_non_inter
