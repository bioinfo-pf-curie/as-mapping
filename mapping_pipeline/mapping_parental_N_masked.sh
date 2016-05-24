#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#       Script to perform alignment to parental genomes
#       STEP 1 : Generation of parental genomes
#       STEP 2 : Masking both parental genomes
#		STEP 3 : Mapping of the reads on both genomes
#       STEP 4 : BAM processing for Allele Specific analysis

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

sam_out=${sam_out}mapping_parental_N_masked/
mkdir -p $sam_out

##### STEP 1 : Generation of parental genomes ------------------------------------------------
#       Following this script, both parental genomes are generated in the $fasta_out directory
#       with the named specified in $fasta_geno1 and $fasta_geno2

if [[ ! -e $fasta_out$fasta_geno1 && ! -e $fasta_out$fasta_geno2 ]]
then
        ${parental_genomes} -c ${config}
fi


##### STEP 2 : Generation of parental N-masked genome ----------------------------------------
#       Following this step, N-masked genome is generated based on differential SNPs between 
#       the two strains specified by $id_geno1 and $îd_geno2

#	Differential SNPs between the two strains
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

#	Generation of masked_genome for both parents (for analysis, should lead to the same masked genome)
geno1_masked_genome="N_masked_"$id_geno1"_"$id_geno2"_from"$id_geno1
geno2_masked_genome="N_masked_"$id_geno1"_"$id_geno2"_from"$id_geno2

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

# Generation of masked genome from geno1
if [[ ! -e $fasta_out$geno1_masked_genome.fa ]]
then # Need to generate the masked genome for geno1
        $genomeMask -r $geno1_masked_genome -f ${fasta_out}chromosomes${id_geno1}/ -v ${diff_bed} -o $fasta_out
        end_genomeMask=`date +%s`
        echo "Masked genome $id_geno1 generated in "$((end_genomeMask-start))" seconds."
fi

# Generation of masked genome from geno2
if [[ ! -e $fasta_out$geno2_masked_genome.fa ]]
then # Need to generate the masked genome for geno2
        $genomeMask -r $geno2_masked_genome -f ${fasta_out}chromosomes${id_geno2}/ -v ${diff_bed} -o $fasta_out
        end_genomeMask=`date +%s`
        echo "Masked genome $id_geno2 generated in "$((end_genomeMask-start))" seconds."
fi


##### STEP 3 : Alignment with Bowtie2 -----------------------------------------------------------

# Updating the location of the bowtie2 indexes
bowtie2_indexes=$fasta_out$bowtie2_indexes

${bowtie2}bowtie2 $SCORING_OPT -p 8 -N 1 -x ${bowtie2_indexes}${geno1_masked_genome} -U ${fq_reads} -S ${sam_out}${geno1_masked_genome}.sam
${bowtie2}bowtie2 $SCORING_OPT -p 8 -N 1 -x ${bowtie2_indexes}${geno2_masked_genome} -U ${fq_reads} -S ${sam_out}${geno2_masked_genome}.sam


# Transform to BAM format and delete SAM file
for masked_genome in ${geno1_masked_genome} ${geno2_masked_genome}
do	
	${samtools} view -bS ${sam_out}${masked_genome}.sam > ${sam_out}${masked_genome}.bam
	rm ${sam_out}${masked_genome}.sam
done

##### STEP 4 : BAM analysis ---------------------------------------------------------------------

mkdir -p ${sam_out}AllelicStatus ${sam_out}mpileup
rm ${sam_out}mpileup/CONFIG
for masked_genome in ${geno1_masked_genome} ${geno2_masked_genome}
do
	#	Script to add a flag with the allelic status for each read
	${markAllelicStatus} -i ${sam_out}${masked_genome}.bam -s ${diff_vcf} -r -o ${sam_out}AllelicStatus/${masked_genome}_withAS.bam
	#   Mpileup to counts the bases present at every SNP positions in the read
	${samtools} sort ${sam_out}${masked_genome}.bam ${sam_out}sorted_${masked_genome} && ${samtools} index ${sam_out}sorted_${masked_genome}.bam
	${samtools} mpileup -l ${diff_bed} -Q 0 ${sam_out}sorted_${masked_genome}.bam > ${sam_out}mpileup/${masked_genome}.pileup
	echo -e "${masked_genome}\t${sam_out}mpileup/${masked_genome}.pileup" >> ${sam_out}mpileup/CONFIG
done
${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_parentalNmask.txt
${sort_counts} -i ${sam_out}mpileup/counts_mapping_parentalNmask.txt -s ${diff_vcf} > ${sam_out}mpileup/counts_parentalNmask.txt

#	Split both counts for each genotype
for masked_genome in ${geno1_masked_genome} ${geno2_masked_genome}
do
	awk -v p=${masked_genome} '{if ($1 ~ p) print}' ${sam_out}mpileup/counts_parentalNmask.txt > ${sam_out}mpileup/final_counts_${masked_genome}.txt
done

#   Cleaning 
rm ${sam_out}sorted_*.b* ${sam_out}mpileup/${masked_genome}.pileup

#   Comparison between generated BAM and mapped reads
gen_bam=${fq_reads%.fq.gz}.bam
if [[ -e ${gen_bam} ]]
then
    mkdir -p ${sam_out}comptoGen_from${id_geno1} ${sam_out}comptoGen_from${id_geno2}
    ${samtools} sort -n ${sam_out}${geno1_masked_genome}.bam ${sam_out}nsorted_${geno1_masked_genome}
	${samtools} sort -n ${sam_out}${geno2_masked_genome}.bam ${sam_out}nsorted_${geno2_masked_genome}
    ${samtools} sort -n ${gen_bam} ${sam_out}nsorted_generated
    ${compMaptoGen} -1 ${id_geno1} -2 ${id_geno2} -g ${sam_out}nsorted_generated.bam -m ${sam_out}nsorted_${geno1_masked_genome}.bam -o ${sam_out}comptoGen_from${id_geno1}/
    ${compMaptoGen} -1 ${id_geno1} -2 ${id_geno2} -g ${sam_out}nsorted_generated.bam -m ${sam_out}nsorted_${geno2_masked_genome}.bam -o ${sam_out}comptoGen_from${id_geno2}/
    # Cleaning
    rm ${sam_out}nsorted_${geno1_masked_genome}.bam ${sam_out}nsorted_${geno2_masked_genome}.bam ${sam_out}nsorted_generated.bam
fi


end=`date +%s`
echo " ============================================================== "
echo "  mapping_parental_N_masked.sh run in "$((end-start))" seconds."
echo " ============================================================== "
