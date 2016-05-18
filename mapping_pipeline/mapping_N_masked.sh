#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  \\\   Script to perform alignment to a N masked genome
#		STEP 1 : Generation of the N-masked genome (masked bases are those that differ between the two studies strains
#		STEP 2 : Mapping of the reads on the N-masked genome
#		STEP 3 : BAM analysis for allele specific studies

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
markAllelicStatus=${map_path}src/Nmask/markAllelicStatus.py

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
sam_out=$sam_out"mapping_N_masked_0_penalty/"
mkdir -p $sam_out

##### STEP 1 : Generation of N-masked genome ----------------------------------------------------

#       Following this script, N-masked genome is generated based on differential SNPs between 
#       the two strains specified by $id_geno1 and $îd_geno2

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
fi

masked_genome="N_masked_"$id_geno1"_"$id_geno2
if [[ ! -e $fasta_out$masked_genome.fa ]]
then # Need to generate the masked genome
	# Transform VCF in bed file
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
    $genomeMask -r $masked_genome -f $ref_dir -v ${diff_bed} -o $fasta_out
	end_genomeMask=`date +%s`
	echo "Masked genome generated in "$((end_genomeMask-start))" seconds."
fi


##### STEP 2 : Alignment with Bowtie2 -----------------------------------------------------------

# Updating the location of the bowtie2 indexes
bowtie2_indexes=$fasta_out$bowtie2_indexes

${bowtie2}bowtie2 $SCORING_OPT --np 0 -p 8 -N 1 -x $bowtie2_indexes$masked_genome -U $fq_reads -S ${sam_out}${masked_genome}.sam

# Transform to BAM format and delete SAM file
${samtools} view -bS ${sam_out}${masked_genome}.sam > ${sam_out}${masked_genome}.bam
rm ${sam_out}${masked_genome}.sam

##### STEP 3 : BAM analysis -----------------------------------------------------------

# Script to add a flag with the allelic status for each read 
mkdir -p ${sam_out}AllelicStatus
${markAllelicStatus} -i ${sam_out}${masked_genome}.bam -s ${diff_vcf} -r -o ${sam_out}AllelicStatus/${masked_genome}_withAS.bam

#   Mpileup to counts the bases present at every SNP positions in the read
mkdir -p ${sam_out}mpileup
${samtools} sort ${sam_out}${masked_genome}.bam ${sam_out}sorted_${masked_genome} && ${samtools} index ${sam_out}sorted_${masked_genome}.bam
${samtools} mpileup -l ${diff_bed} -Q 0 ${sam_out}sorted_${masked_genome}.bam > ${sam_out}mpileup/${masked_genome}.pileup
echo -e "mapping_N_masked\t${sam_out}mpileup/${masked_genome}.pileup" > ${sam_out}mpileup/CONFIG
${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_Nmask.txt
${sort_counts} -i ${sam_out}mpileup/counts_mapping_Nmask.txt -s ${diff_vcf} > ${sam_out}mpileup/final_counts_mapping_Nmask.txt

#   Cleaning 
rm ${sam_out}sorted_${masked_genome}.b* ${sam_out}mpileup/${masked_genome}.pileup

end=`date +%s`
echo " ======================================================= "
echo "  mapping_N_masked.sh run in "$((end-start))" seconds."
echo " ======================================================= "
