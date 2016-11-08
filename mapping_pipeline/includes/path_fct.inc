
############################## LOCAL SCRIPTS ##################################

# Building genomes and indexes 
BUILDING_PATH=${PIPELINE_PATH}/scripts/build_genomes
BUILD_B2_INDEX=${BUILDING_PATH}/build_bowtie2_index.sh
BUILD_NMASK_GENOME=${BUILDING_PATH}/build_N-masked_genome.sh
BUILD_PARENTAL_GENOME=${BUILDING_PATH}/build_parental_genome.sh
BUILD_DIPLOID_GENOME=${BUILDING_PATH}/build_diploid_genome.sh

# Mapping scripts
MAPPING_PATH=${PIPELINE_PATH}/scripts/mapping
# - Reference mapping
MAPPING_REFERENCE_B2=${MAPPING_PATH}/bowtie2_mapping_reference.sh
MAPPING_REFERENCE_TOPHAT=${MAPPING_PATH}/tophat_mapping_reference.sh
# - N-masked mapping
MAPPING_NMASKED_B2=${MAPPING_PATH}/bowtie2_mapping_N-masked.sh
MAPPING_NMASKED_TOPHAT=${MAPPING_PATH}/tophat_mapping_N-masked.sh
# - Parental mapping
MAPPING_PARENTAL_B2=${MAPPING_PATH}/bowtie2_mapping_parental.sh
MAPPING_PARENTAL_TOPHAT=${MAPPING_PATH}/tophat_mapping_parental.sh
# - Diploid mapping
MAPPING_DIPLOID_B2=${MAPPING_PATH}/bowtie2_mapping_diploid.sh
MAPPING_DIPLOID_TOPHAT=${MAPPING_PATH}/tophat_mapping_diploid.sh

# Scripts
SCRIPTS_PATH=${PIPELINE_PATH}/scripts
MERGE_ALIGN=${SCRIPTS_PATH}/mergeAlign.py
MARKALLELICSTATUS=${SCRIPTS_PATH}/markAllelicStatus.py

# Not used anymore
EXTRACT_SNPS=${SCRIPTS_PATH}/extracts_snps.py
SNPSPLIT_GENOMES=${SCRIPTS_PATH}/SNPsplit_genomes.sh
DIPLOID_GENOMES=${SCRIPTS_PATH}/diploid_genomes.sh
ANNOTATE_COUNTS=${SCRIPTS_PATH}/annotate_counts.py
COMPMAPTOGEN=${SCRIPTS_PATH}/compMaptoGen.py
SELECT_FROM_DIP=${SCRIPTS_PATH}/dip/select_from_diploid.py

# SNPsplit
SNPSPLIT=${SNPSPLIT_DIR}/SNPsplit
SNPSPLIT_GEN=${SNPSPLIT_DIR}/SNPsplit_genome_preparation

# Samtools
SAMTOOLS=${SAMTOOLS_DIR}/samtools

# Folders for different file type
FASTA_OUT=${GEN_OUT}/fastas
INDEXES=${FASTA_OUT}/indexes
VCF_OUT=${GEN_OUT}/vcfs