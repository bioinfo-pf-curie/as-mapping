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
sam_out=$sam_out"mapping_parental/"

# Select best alignment and mark Allelic status
#   Sort bam files by names
${samtools} sort -n ${sam_out}${id_geno1}.bam ${sam_out}nsorted_${id_geno1}
${samtools} sort -n ${sam_out}${id_geno2}.bam ${sam_out}nsorted_${id_geno2}
${selectBest} -1 ${sam_out}nsorted_${id_geno1}.bam -2 ${sam_out}nsorted_${id_geno2}.bam -f selected_${id_geno1}_${id_geno2} -o ${sam_out}
#	Cleaning files
rm ${sam_out}nsorted_${id_geno1}.bam ${sam_out}nsorted_${id_geno2}.bam

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
echo -e "par\t${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup" > ${sam_out}mpileup/CONFIG

${samtools} mpileup -l ${diff_bed} -Q 0 -q 20 ${sam_out}sorted_${id_geno1}_${id_geno2}.bam > ${sam_out}mpileup/${id_geno1}_${id_geno2}_filtered.pileup
echo -e "par_filtered\t${sam_out}mpileup/${id_geno1}_${id_geno2}_filtered.pileup" >> ${sam_out}mpileup/CONFIG

diff_bed_ori=/bioinfo/users/khillion/AS_proj/generated_reads/0704_adapt_SNPsplit/data/vcfs/diff_SNPs.bed
${samtools} mpileup -l ${diff_bed_ori} -Q 0 ${sam_out}sorted_${id_geno1}_${id_geno2}.bam > ${sam_out}mpileup/${id_geno1}_${id_geno2}_origin.pileup
echo -e "par_origin\t${sam_out}mpileup/${id_geno1}_${id_geno2}_origin.pileup" >> ${sam_out}mpileup/CONFIG

${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_parental.txt
${annotate_counts} -i ${sam_out}mpileup/counts_mapping_parental.txt -s ${diff_vcf} > ${sam_out}mpileup/final_counts_mapping_parental.txt
#   Cleaning 
rm ${sam_out}sorted_${id_geno1}_${id_geno2}.b* ${sam_out}non_ambiguous.bam ${sam_out}mpileup/*.pileup

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
