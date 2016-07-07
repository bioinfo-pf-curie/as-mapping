#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :

#### Parameters #### --------------------------------------------------------------------

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

# -- Local scripts for BAM analysis
selectBest=${map_path}scripts/par/select_best.py

#### Function #### ----------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------

# Set up output directory
sam_out=${sam_out}mapping_diploid/

# First step is to rename the header and the chromosome names
${samtools} view -H ${sam_out}${id_geno1}_${id_geno2}.bam | grep -v -E "chr[A-Z0-9]+_${id_geno2}" | sed "s/_${id_geno1}\s/\t/" > ${sam_out}SAM.tmp
${samtools} view ${sam_out}${id_geno1}_${id_geno2}.bam | sed -E 's/_[0-9a-Z_]+\s/\t/' >> ${sam_out}SAM.tmp
${samtools} view -bS ${sam_out}SAM.tmp > ${sam_out}${id_geno1}_${id_geno2}_renamed.bam
rm ${sam_out}SAM.tmp 

${select_from_dip} -b ${sam_out}${id_geno1}_${id_geno2}_renamed.bam -n selected_${id_geno1}_${id_geno2} -o ${sam_out} 

rm ${sam_out}${id_geno1}_${id_geno2}_renamed.bam

# Build different VCF and BED files for analysis
mkdir -p ${vcf_out}
vcf=$(basename $full_vcf)
diff_vcf=$vcf_out${vcf%.vcf}"_"$id_geno1"_"$id_geno2".vcf"
if [[ ! -e $diff_vcf ]]
then # Need to generate the VCF with the differential SNPs
    echo "Generating VCF file of different SNPs between $id_geno1 and $id_geno2 ..."
    $extract_SNPs -i $full_vcf -r $id_geno1 -a $id_geno2 -f 1 > $diff_vcf
fi
diff_bed=$vcf_out${vcf%.vcf}"_"$id_geno1"_"$id_geno2".bed"
if [[ ! -e $diff_bed ]]
then # Need to generate the BED with the differential SNPs from VCF
    echo "Generating BED file of different SNPs between $id_geno1 and $id_geno2 from VCF file ..."
    awk '{if($1 !~ /^#/) print "chr"$1,$2-1,$2,$3":"$4"/"$5}' OFS='\t' $diff_vcf > $diff_bed
fi

# Mpileup to counts the bases present at every SNP positions in the read
mkdir -p ${sam_out}mpileup
${samtools} view -h ${sam_out}selected_${id_geno1}_${id_geno2}.bam | grep -v 'XA:i:3' | ${samtools} view -bS - > ${sam_out}non_ambiguous.bam
${samtools} sort ${sam_out}non_ambiguous.bam ${sam_out}sorted_${id_geno1}_${id_geno2} && ${samtools} index ${sam_out}sorted_${id_geno1}_${id_geno2}.bam
${samtools} mpileup -l ${diff_bed} -Q 0 ${sam_out}sorted_${id_geno1}_${id_geno2}.bam > ${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup
echo -e "mapping_diploid\t${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup" > ${sam_out}mpileup/CONFIG
${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_diploid.txt
${annotate_counts} -i ${sam_out}mpileup/counts_mapping_diploid.txt -s ${diff_vcf} > ${sam_out}mpileup/final_counts_mapping_diploid.txt
#   Cleaning 
rm ${sam_out}sorted_${id_geno1}_${id_geno2}.b* ${sam_out}non_ambiguous.bam ${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup

# Comparison between generated BAM and mapped reads
gen_bam=${fq_reads%.fq.gz}.bam
if [[ -e ${gen_bam} ]]
then
    mkdir -p ${sam_out}comptoGen
    ${samtools} sort -n ${sam_out}selected_${id_geno1}_${id_geno2}.bam ${sam_out}nsorted_selected
    ${samtools} sort -n ${gen_bam} ${sam_out}nsorted_generated
    ${compMaptoGen} -1 ${id_geno1} -2 ${id_geno2} -g ${sam_out}nsorted_generated.bam -m ${sam_out}nsorted_selected.bam -o ${sam_out}comptoGen/ -u
    # Cleaning
    rm ${sam_out}nsorted_selected.bam ${sam_out}nsorted_generated.bam
fi
