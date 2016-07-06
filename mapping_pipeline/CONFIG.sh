#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) : CONFIG file for mapping of reads with different methods :
#	-1 : Mapping to reference genome
#	-2 : Mapping to parental genomes (both chromosomes separatly)
#	-3 : Mapping to diploid genome (both parental chromosomes in the same time)
#	-4 : Mapping to reference N-masked genome
#	-5 : Mapping to both parental N-masked genomes

############################## CONFIG FILE ##################################


# ------------ Tools path(s) ------------

# Required installed tools
bowtie2=/bioinfo/local/build/bowtie2/bowtie2-2.2.5/
vcf2diploid=/bioinfo/users/khillion/bin/vcf2diploid_v0.2.6a/vcf2diploid.jar # Obsolete with SNPsplit
samtools=/bioinfo/local/build/samtools/samtools/bin/samtools
checkVariants=/bioinfo/users/khillion/bin/clinTools/check_variants.py
SNPsplit_gen=/bioinfo/users/khillion/GIT/SNPsplit/SNPsplit_genome_preparation
# Local scripts
#	path to the directory with the different scripts
map_path="/bioinfo/users/khillion/GIT/as-mapping/mapping_pipeline/"
extract_SNPs=${map_path}scripts/extract_snps.py
genomeMask=${map_path}scripts/genomeMask.sh             # Obsolete with SNPsplit
parental_genomes=${map_path}scripts/parental_genomes.sh # Obsolete with SNPsplit
SNPsplit_genomes=${map_path}scripts/SNPsplit_genomes.sh
diploid_genome=${map_path}scripts/diploid_genome.sh
annotate_counts=${map_path}scripts/annotate_counts.py
compMaptoGen=${map_path}scripts/compMaptoGen.py
select_from_dip=${map_path}scripts/dip/select_from_diploid.py

 
# ------------ Input(s) ------------

# -- Required for all Mappings
fq_reads="READS_PATH"
ref_geno="/data/annotations/Mouse/mm10/complete/mm10.fa"

# -- Required for Mapping 2 to 5
# 	You can specify "ref" for the first genome
id_geno1="CAST_EiJ"
id_geno2="129S1_SvImJ"

# -- Required for Mapping 2,3 and 5
#	VCF containing strain specific SNPs
#	Can be downloaded at for v5 : ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs/
vcf_geno1="/data/annotations/Mouse/variant_informations/mgpV5_mm10/strain_specific_vcfs_SANGER/CAST_EiJ/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
vcf_geno2="/data/annotations/Mouse/variant_informations/mgpV5_mm10/strain_specific_vcfs_SANGER/129S1_SvImJ/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf"

# -- Required for Mapping 4 and 5
#	Directory containing every chromosomes of the reference genome at the fasta format
ref_dir="/data/annotations/Mouse/mm10/chromosomes_clean/"

#       full_vcf : vcf file containing all SNPs from either the two strains, 
#		or vcf file from Sanger directly (all SNPs of all strains, takes more time to extract information)
full_vcf="/data/annotations/Mouse/variant_informations/mgpV5_mm10/mgp.v5.merged.snps_all.dbSNP142.vcf"
#diff_vcf="/bioinfo/users/khillion/AS_proj/mapping_pipeline/data/vcfs/mgp.v5.merged.snps_all.dbSNP142_CAST_EiJ_129S1_SvImJ.vcf"
diff_bed="/bioinfo/users/khillion/AS_proj/mapping_pipeline/data/vcfs/mgp.v5.merged.snps_all.dbSNP142_CAST_EiJ_129S1_SvImJ.bed"



# ------------ Output(s) ------------

# -- Directories

# 	Main output directories for files to be reused (N masked genome, vcf files ...)
main_out="/bioinfo/users/khillion/AS_proj/mapping_pipeline/data/"
#	Folders for different file type
fasta_out=${main_out}"fastas/"
vcf_out=${main_out}"vcfs/"
#	Name of the folder containing bowtie2 indexes, MUST BE in the same folder as the .fa file (only specified folder name here)
bowtie2_indexes="bowtie2_indexes/"
#	Output directory for BAM
sam_out="MAPPED_OUTDIR"


# -- Files (Do not leave empty)

#	Names of VCF files containing only homozygous SNPs of a strain
#homo_vcf1="homo_"$(echo $vcf_geno1 | tr '/' '\n' | tail -1)
#homo_vcf2="homo_"$(echo $vcf_geno2 | tr '/' '\n' | tail -1)
#	Names of fasta files of parental genomes
#fasta_geno1=$id_geno1".fa"
#fasta_geno2=$id_geno2".fa"


# ------------ Bowtie2 Options ------------

# -- Scoring options (Current are default).
#	--np 1			# --np : penalty for matching to an ambiguous character such as "N"
#	--mp 6,2		# --mp MX,MN : maximum and minimum mismatch penalty
#	--rdg 5,3		# --rdg : sets the read gap open and extend
#	--rfg 5,3		# --rfg : sets the reference gap open and extend
#	--score-min L,-0.6,-0.6	# --score-min : function for the minimum alignement score to be considered valid
SCORING_OPT=""
