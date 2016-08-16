
############################## CONFIG FILE ##################################

# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : CONFIG file for mapping of reads with different methods :
#	 : Mapping to reference genome
#	 : Mapping to parental genomes (both chromosomes separatly)
#	 : Mapping to diploid genome (both parental chromosomes in the same time)
#	 : Mapping to reference N-masked genome

# ------------ Tool paths ------------

# Required installed tools
PYTHON=/bioinfo/local/build/python/python-2.7.9/bin/python2.7
BOWTIE2_DIR=/bioinfo/local/build/bowtie2/bowtie2-2.2.5/
SAMTOOLS=/bioinfo/local/build/samtools/samtools/bin/samtools
CHECK_VARIANTS=/bioinfo/users/khillion/bin/clinTools/check_variants.py
SNPSPLIT_GEN=/bioinfo/users/khillion/GIT/SNPsplit/SNPsplit_genome_preparation


# Local scripts
#	Paths to the directory with the different scripts
MAP_PATH=/bioinfo/users/khillion/GIT_dev/as-mapping/mapping_pipeline
#   Get local scripts from the pipeline
source LOCAL_SCRIPTS ${MAP_PATH}

 
# ------------ Input(s) ------------

# Files or Directories
#   Path of the reads to map (either .fq or .fq.gz)
FQ_READS="READS_PATH"
#   Full reference genome (fasta format)
ref_geno="/data/annotations/Mouse/mm10/complete/mm10.fa"
#   Directory containing every chromosomes of the reference genome at the fasta format
#   
ref_dir="/data/annotations/Mouse/mm10/chromosomes_clean/"
#   full_vcf : vcf file containing all SNPs from either the two strains, 
full_vcf="/data/annotations/Mouse/variant_informations/mgpV5_mm10/mgp.v5.merged.snps_all.dbSNP142.vcf"

# Informations
# 	Genotype names
id_geno1="CAST_EiJ"
id_geno2="129S1_SvImJ"


# ------------ Output(s) ------------

# Output directories for files to be reused (N masked genome, indexes ...)
main_out="/bioinfo/users/khillion/AS_proj/mapped_reads/data/"

#	Output directory for BAM
sam_out="MAPPED_OUTDIR"

# OTHER CONFIG FILE ?
#	Folders for different file type
fasta_out=${main_out}"fastas/"
vcf_out=${main_out}"vcfs/"
#	Name of the folder containing bowtie2 indexes, MUST BE in the same folder as the .fa file (only specified folder name here)
bowtie2_indexes="${fasta_out}bowtie2_indexes/"


# ------------ Bowtie2 Options ------------

# -- Main options (non exhaustive list)
#	-p 1   			#
#	--reorder		#
#	--sensitive		#
B2_OPTIONS="--reorder -p 8 --sensitive"

# -- Scoring options (Current are default).
#	--np 1			# --np : penalty for matching to an ambiguous character such as "N"
#	--mp 6,2		# --mp MX,MN : maximum and minimum mismatch penalty
#	--rdg 5,3		# --rdg : sets the read gap open and extend
#	--rfg 5,3		# --rfg : sets the reference gap open and extend
#	--score-min L,-0.6,-0.6	# --score-min : function for the minimum alignement score to be considered valid
B2_SCORING_OPT=""
