#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#  \\\	Script to map reads on reference genome using Bowtie2
#		STEP 1 : Mapping of the reads
#		STEP 2 : BAM analysis for Allele specific analysis


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
markAllelicStatus=${map_path}src/ref/markAllelicStatus.py

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

##### STEP 1 : Alignment to reference genome ----------------------------------------------------

# Set up output directory for this method of mapping in the main output directory
sam_out=${sam_out}mapping_reference/
mkdir -p ${sam_out}

# Get the name of the reference genome to use already existing indexes (might change this part to work directly on my own indexes)
id_ref=$(basename ${ref_geno})
id_ref=${id_ref%.fa}

# Updating the location of the bowtie2 indexes (might change as well, same reason as above)
bowtie2_indexes=${ref_geno%$id_ref.fa}${bowtie2_indexes}

# Create bowtie2 indexes if they do not exist
if [ ! -e ${bowtie2_indexes}${id_ref}.rev.2.bt2 ]
then
	echo "Generating bowtie2 indexes ..."
	mkdir $bowtie2_indexes
	${bowtie2}bowtie2-build -f ${fasta_out}${id_ref}.fa ${bowtie2_indexes}${id_ref}
fi

#${bowtie2}bowtie2 $SCORING_OPT -p 8 -N 1 -x ${bowtie2_indexes}${id_ref} -U $fq_reads -S ${sam_out}${id_ref}".sam"

# Transform to BAM format and delete SAM file
#${samtools} view -bS ${sam_out}${id_ref}.sam > ${sam_out}${id_ref}.bam
#rm ${sam_out}${id_ref}.sam

##### STEP 2 : BAM analysis ----------------------------------------------------

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

# 	Script to add a flag with the allelic status for each read
 
mkdir -p ${sam_out}AllelicStatus
#${markAllelicStatus} -i ${sam_out}${id_ref}.bam -s ${diff_vcf} -r -o ${sam_out}AllelicStatus/${id_ref}_withAS.bam

#	Mpileup to counts the bases present at every SNP positions in the read
mkdir -p ${sam_out}mpileup
#${samtools} sort ${sam_out}${id_ref}.bam ${sam_out}sorted_${id_ref} && ${samtools} index ${sam_out}sorted_${id_ref}.bam
#${samtools} mpileup -l ${diff_bed} -Q 0 ${sam_out}sorted_${id_ref}.bam > ${sam_out}mpileup/${id_ref}.pileup
echo -e "mapping_reference\t${sam_out}mpileup/${id_ref}.pileup" > ${sam_out}mpileup/CONFIG
#${checkVariants} ${sam_out}mpileup/CONFIG ${ref_geno} > ${sam_out}mpileup/counts_mapping_reference.txt
#${sort_counts} -i ${sam_out}mpileup/counts_mapping_reference.txt -s ${diff_vcf} > ${sam_out}mpileup/final_counts_mapping_reference.txt
# 	Cleaning 
#rm ${sam_out}sorted_${id_ref}.b* ${sam_out}mpileup/${id_ref}.pileup

#	Comparison between generated BAM and mapped reads
gen_bam=${fq_reads%.fq.gz}.bam
if [[ -e ${gen_bam} ]]
then
	mkdir -p ${sam_out}comptoGen
	${samtools} sort -n ${sam_out}${id_ref}.bam ${sam_out}nsorted_${id_ref}
	${samtools} sort -n ${gen_bam} ${sam_out}nsorted_generated
	${compMaptoGen} -1 ${id_geno1} -2 ${id_geno2} -g ${sam_out}nsorted_generated.bam -m ${sam_out}nsorted_${id_ref}.bam -o ${sam_out}comptoGen/
	# Cleaning
	rm ${sam_out}nsorted_${id_ref}.bam ${sam_out}nsorted_generated.bam
fi

end=`date +%s`
echo " ======================================================= "
echo "	mapping_reference.sh run in "$((end-start))" seconds."
echo " ======================================================= "
