#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
# Script to generate bed files from outputs of ART and rephase them to the reference genome

#### Parameters ####
source config.sh

#### Function ####


#### Main ####
start=`date +%s`

# 1) Generate bed files from .aln output of ART
echo "aln2bed is running ..."

geno1_aln=`ls $art_outdir$id_geno1*.aln | tr '/' '\n' | tail -1`
geno2_aln=`ls $art_outdir$id_geno2*.aln | tr '/' '\n' | tail -1`
art_bed_geno1=${geno1_aln%.aln}".bed"
art_bed_geno2=${geno2_aln%.aln}".bed"

#$aln2bed $bed_outdir$art_bed_geno1 $art_outdir$geno1_aln 
#$aln2bed $bed_outdir$art_bed_geno2 $art_outdir$geno2_aln

aln2bed_time=`date +%s`
echo "aln2bed run in "$((aln2bed_time-start))" seconds"


# 2) Rephasing bed files from ART's generated reads using awk
echo "Rephasing of the bed file with awk ..."

phased_geno1="${art_bed_geno1%.fa.bed}.phased.bed"
phased_geno2="${art_bed_geno2%.fa.bed}.phased.bed"

#awk '{split($1,a,":"); split(a[2],b,"-"); print a[1],$2+b[1],$3+b[1],$4,$5,$6}' OFS='\t' $bed_outdir$art_bed_geno1 > $bed_outdir$phased_geno1
#awk '{split($1,a,":"); split(a[2],b,"-"); print a[1],$2+b[1],$3+b[1],$4,$5,$6}' OFS='\t' $bed_outdir$art_bed_geno2 > $bed_outdir$phased_geno2

awk_time=`date +%s`
echo "Awk run in "$((awk_time-start))" seconds"


# 3) Output coverage for both simulated set

# a. Generation of the bed file from the vcf file containing differences between two genotypes
echo "Generating bed file from VCF ..."

if [ -z $diff_vcf ]
then
	vcf=`echo $full_vcf | tr '/' '\n' | tail -1`
        diff_vcf=$vcf_outdir${vcf%.vcf}$id_geno1"_"$id_geno2".vcf"	
fi
vcf_bed=`echo $diff_vcf | tr '/' '\n' | tail -1 | sed s/\.vcf/\.bed/`
awk '{if ($1 !~ /^#/) print $1,$2-1,$2,$3":"$4"/"$5}' OFS="\t" $diff_vcf > $bed_outdir$vcf_bed

# b. Generation of bedfile coverage for both genomes
echo "Bedtools coverage is running ..."

cov_ref="${phased_geno1%.bed}_coverage.bed"
cov_alt="${phased_geno2%.bed}_coverage.bed"

$bedtools coverage -a $bed_outdir$phased_geno1 -b $bed_outdir$vcf_bed -d > $bed_outdir$cov_ref
$bedtools coverage -a $bed_outdir$phased_geno2 -b $bed_outdir$vcf_bed -d > $bed_outdir$cov_alt

bedtools_time=`date +%s`
echo "Bedtools coverage run in "$((bedtools_time-start))" seconds"

end=`date +%s`
echo " ==================================== "
echo "Script run in "$((end-start))" seconds."
echo " ==================================== "
