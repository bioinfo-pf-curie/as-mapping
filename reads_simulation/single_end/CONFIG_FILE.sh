#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :
#	Config file for reads simulation using ART

############################## CONFIG FILE ##################################

# ------------ Path(s) for the tools ------------

# 	Path of the required tools
bedtools=/bioinfo/local/build/BEDTools/BEDTools_2.21.0/bin/bedtools
art=/bioinfo/users/khillion/Tools/ART/art/art_illumina
samtools=/bioinfo/local/build/samtools/samtools/bin/samtools
SNPsplit_gen=/bioinfo/users/khillion/GIT/SNPsplit/SNPsplit_genome_preparation

# 	Path of the read simulation scripts
simreads=/bioinfo/users/khillion/GIT/as-mapping/reads_simulation/single_end/
extract_SNPs=${simreads}scripts/extract_snps.py


# ------------ INPUT ------------

#	IDs of the strains
id_geno1="CAST_EiJ"
id_geno2="129S1_SvImJ"

#   Coverage
coverage=100

#   BED file with the desired region for generation
regions="REGIONS_SELECTED.bed"

#   Directory containing every chromosomes of the reference genome at the fasta format (no "chr")
ref_dir="/bioinfo/users/khillion/data/genomes/mm10/chromosomes/"

# 	VCF files
#		 full_vcf : Contain SNPs of all different strains of mice
full_vcf="/data/annotations/Mouse/variant_informations/mgpV5_mm10/mgp.v5.merged.snps_all.dbSNP142.vcf"

#	[OPTIONAL] Mappability track :
#	Mappability tracks can be generated using GEMTOOLS (gemtools-1.7.1-i3)
#	available at : https://github.com/gemtools/gemtools
mappa="/data/annotations/Mouse/mappability/mm10/100bp/mm10_100.bed"

#	[OPTIONAL] BED file for allele ratio of genotype 1
#	BED5 file with the following field
#	Exemple :
#	CHR	begin	end	name	ratio_geno1
#	chrX	103460372	103483233	Xist	1
ASratio=""


# ------------ ART parameters ------------

read_length=100
subs=0	 			# 93 for no substitution, 0 for default rate of substitution, -xx to increase substitution
rs=1				# Random seed, it is fixed to 1 to use the same for both generation.
ir=0				# Insertion rate (default: 0.00009)
dr=0				# Deletion rate (default: 0.00011)
sequencer="HS20"	# The name of Illumina sequencing system :
					#			GA1  - Genome Analyzer I
                	#			GA2  - Genome Analyzer II
                	#			HS10 - HiSeq 1000
                	#			HS20 - HiSeq 2000
                	#			HS25 - HiSeq 2500 
                	#			MS   - MiSeq


# ------------ OUTPUT ------------

#	Directories 
#		Output directory where new generated files that can be re-used are stored
main_out="data/"
vcf_outdir=${main_out}"vcfs/"
fasta_outdir=${main_out}"fastas/"
tmp_outdir=${main_out}"tmp/"

#		Directory where the generated reads will be stored
art_outdir="simulated_reads/"
