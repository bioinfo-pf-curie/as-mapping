#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  \\\	 Script to perform alignment to parental genomes
#		STEP 1 : Generation of parental genomes
#		STEP 2 : Mapping of the reads on both genomes
#		STEP 3 : Select best alignement between both genomes and BAM analysis for Allele specific analysis

#### Parameters #### --------------------------------------------------------------------------------------------------------------

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
selectBest=${map_path}src/par/select_best.py

#### Function #### ----------------------------------------------------------------------------------------------------------------

# Get args
function usage {
    echo -e "Usage : $0"
    echo -e "-c"" <Config file>"
    echo -e "-h"" <help>"
    exit
}

#### Main #### --------------------------------------------------------------------------------------------------------------------

start=`date +%s`

# Set up output directory for this method of mapping in the main output directory
sam_out=$sam_out"mapping_parental/"
mkdir -p $sam_out


##### STEP 1 : Generation of parental genomes ------------------------------------------------
#	Following this script, both parental genomes are generated in the $fasta_out directory
#	with the named specified in $fasta_geno1 and $fasta_geno2

if [[ ! -e $fasta_out$fasta_geno1 && ! -e $fasta_out$fasta_geno2 ]]
then # Create parental genomes and indexes (generation of indexes within the script)
	${parental_genomes} -c ${config}
fi


##### STEP 2 : Alignment with Bowtie2 --------------------------------------------------------

# Updating the location of the bowtie2 indexes
bowtie2_indexes=$fasta_out$bowtie2_indexes

${bowtie2}bowtie2 $SCORING_OPT -p 8 -N 1 -x $bowtie2_indexes$id_geno1 -U $fq_reads -S ${sam_out}${id_geno1}.sam
${bowtie2}bowtie2 $SCORING_OPT -p 8 -N 1 -x $bowtie2_indexes$id_geno2 -U $fq_reads -S ${sam_out}${id_geno2}.sam

# Transform to BAM format and delete SAM file
${samtools} view -bS ${sam_out}${id_geno1}.sam > ${sam_out}${id_geno1}.bam && rm ${sam_out}${id_geno1}.sam
${samtools} view -bS ${sam_out}${id_geno2}.sam > ${sam_out}${id_geno2}.bam && rm ${sam_out}${id_geno2}.sam


##### STEP 3 : Comparison to find best alignment and BAM analysis  ---------------------------

# Select best alignment and mark Allelic status
#	Sort bam files by names
${samtools} sort -n ${sam_out}${id_geno1}.bam ${sam_out}nsorted_${id_geno1}
${samtools} sort -n ${sam_out}${id_geno2}.bam ${sam_out}nsorted_${id_geno2}
${selectBest} -1 ${sam_out}nsorted_${id_geno1}.bam -2 ${sam_out}nsorted_${id_geno2}.bam -f selected_${id_geno1}_${id_geno2} -o ${sam_out}
rm ${sam_out}nsorted_${id_geno1}.bam ${sam_out}nsorted_${id_geno2}.bam

if [[ -z $diff_vcf ]]
then # User did not specify a VCF with differential SNPs
    mkdir -p ${vcf_out}
    vcf=$(basename $full_vcf)
    diff_vcf=$vcf_out${vcf%.vcf}"_"$id_geno1"_"$id_geno2".vcf"
    if [[ ! -e $diff_vcf ]]
    then # Need to generate the VCF with the differential SNPs
        echo "Generating VCF file of different SNPs between $id_geno1 and $id_geno2 ..."
        $extract_SNPs -i $full_vcf -r $id_geno1 -a $id_geno2 -f 1 > $diff_vcf
        end_vcf_gen=`date +%s`
        echo "VCF file generated in "$((end_vcf_gen-start))" seconds."
    fi  
    if [[ -z $diff_bed ]]
    then # User did not specify a BED with differential SNPs
        diff_bed=$vcf_out${vcf%.vcf}"_"$id_geno1"_"$id_geno2".bed"
        if [[ ! -e $diff_bed ]]
        then # Need to generate the BED with the differential SNPs from VCF
            echo "Generating BED file of different SNPs between $id_geno1 and $id_geno2 from VCF file ..."
            awk '{if($1 !~ /^#/) print "chr"$1,$2-1,$2,$3":"$4"/"$5}' OFS='\t' $diff_vcf > $diff_bed
            end_bed_gen=`date +%s`
            echo "BED file generated in "$((end_bed_gen-start))" seconds."
        fi  
    fi  
fi

#   Mpileup to counts the bases present at every SNP positions in the read
mkdir -p ${sam_out}mpileup
${samtools} sort ${sam_out}selected_${id_geno1}_${id_geno2}.bam ${sam_out}sorted_${id_geno1}_${id_geno2} && ${samtools} index ${sam_out}sorted_${id_geno1}_${id_geno2}.bam
${samtools} mpileup -l ${diff_bed} -Q 0 ${sam_out}sorted_${id_geno1}_${id_geno2}.bam > ${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup
echo -e "mapping_parental\t${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup" > ${sam_out}mpileup/CONFIG
${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_parental.txt
${sort_counts} -i ${sam_out}mpileup/counts_mapping_parental.txt -s ${diff_vcf} > ${sam_out}mpileup/final_counts_mapping_parental.txt

#   Cleaning 
rm ${sam_out}sorted_${id_geno1}_${id_geno2}.b* ${sam_out}mpileup/${id_geno1}_${id_geno2}.pileup

end=`date +%s`
echo " ====================================================== "
echo "  mapping_parental.sh run in "$((end-start))" seconds."
echo " ====================================================== "
