#!/bin/bash
# Author(s) : Kenzo Hillion
# Contact : kenzo.hillion@curie.fr
# Comment(s) :

############################## CONFIG FILE ##################################

# ------------ Path(s) for the tools ------------
bedtools=/bioinfo/local/build/BEDTools/BEDTools_2.21.0/bin/bedtools
art=/bioinfo/users/khillion/Tools/ART/art/art_illumina				# Configuration for ART at the end of this config file
vcf2diploid=/bioinfo/users/khillion/bin/vcf2diploid_v0.2.6a/vcf2diploid.jar
extract_SNPs=/bioinfo/users/khillion/bin/extract_snps_forNmask.py
samtools=/bioinfo/local/build/samtools/samtools/bin/samtools
simreads=/bioinfo/users/khillion/GIT/as-mapping/reads_simulation/single_end/


# ------------ INPUT ------------

#	IDs of the strains
id_geno1="CAST_EiJ"
id_geno2="129S1_SvImJ"

# 	Fasta of the reference chromosome
ref="/data/annotations/Mouse/mm10/chromosomes/chrX.fa"
#	\__	Name of the working chromosome
chr=$(basename ${ref%.fa}) && chr=${chr:3}

# 	VCF files
#		 diff_vcf : Contain SNPs which are different between the two strains (will be generated if left empty)
#		 full_vcf : Contain SNPs of all different strains of mice
vcf_geno1="/data/annotations/Mouse/variant_informations/mgpV5_mm10/strain_specific_vcfs_SANGER/CAST_EiJ/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
vcf_geno2="/data/annotations/Mouse/variant_informations/mgpV5_mm10/strain_specific_vcfs_SANGER/129S1_SvImJ/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf"
diff_vcf=""
full_vcf="/data/annotations/Mouse/variant_informations/mgpV5_mm10/chromosomes/chrX_mgp.v5.merged.snps_all.dbSNP142.vcf"

#	Mappability track (optional) :
mappa="/data/annotations/Mouse/mappability/mm10/100bp/mm10_100.bed"

#	BED files of chosen genes/regions for allele specificity :
ASspecgeno1=""
ASspecgeno2=""


# ------------ OUTPUT ------------

#	Directories 
#		Output directory where new generated files that can be re-used are stored
main_out="data/"
vcf_outdir=${main_out}"vcfs/"
fasta_outdir=${main_out}"fastas/"
bed_outdir=${main_out}"beds/"

#		Directory where the generated reads will be stored
art_outdir="OUTPUT_NAME/"

# 	BED file of selected intervals for read generation (which is reused in case of several simulation)
int_bed="intervals_100.bed"


# ------------ ART parameters ------------

read_length=100
subs=93	 			# 93 for no substitution, 0 for default rate of substitution (seems to be 0.0074/base) -xx to increase substitution
rs=1				# Random seed, it is fixed to 1 to use the same for both generation.
ir=0				# Insertion rate
dr=0				# Deletion rate
sequencer="HS20"	# The name of Illumina sequencing system :
					#			GA1  - Genome Analyzer I
                	#			GA2  - Genome Analyzer II
                	#			HS10 - HiSeq 1000
                	#			HS20 - HiSeq 2000
                	#			HS25 - HiSeq 2500 
                	#			MS   - MiSeq
