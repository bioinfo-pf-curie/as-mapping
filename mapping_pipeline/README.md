# Mapping Pipeline for Allele Specific analysis in the Mouse

------------------------------

This page is a quick start guide, please read the full [online manual](link) for more information.

See NEWS for information about changes in this and previous versions.

## What is it ?

This pipeline was designed for the mapping of reads prior to allele specific analysis. Most of the described methods are present in the present pipeline (reference genome, parental genomes, diploid genome and N-masked genome).

## Contact

For any question about the pipeline, please contact <kenzo.hillion@curie.fr>.

------------------------------

# Quick Start Guide

## Requirements

Python version 2.7 is used for this pipeline and the [Pysam](https://github.com/pysam-developers/pysam) (version 0.8.4) Python module is required.

#### Dependencies:

* [SNPsplit](http://github.com/FelixKrueger/SNPsplit) (version 0.3)
* [samtools](http://samtools.sourceforge.net) (version 1.1)

#### Mappers:

The current version of the pipeline supports mapping with either Bowtie2 or Tophat:

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 2.2.5)
* [Tophat](http://ccb.jhu.edu/software/tophat/index.shtml) (version 2.1.0)


## CONFIG File

The `CONFIG` file contains all the informations and paths for a simulation.

#### The different path of the tools:

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
PYTHON            | path to Python
SAMTOOLS          | path to samtools
SNPSPLIT_DIR      | path to SNPsplit directory
PIPELINE_PATH     | path of the mapping pipeline


### Inputs:

#### Reference files and informations:

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
FULL_VCF          | path to the VCF file (.vcf). For mouse, VCF files are downloaded [here](ftp://ftp-mouse.sanger.ac.uk/)
REF_GENO          | reference genome (.fa)
B2_INDEX_REF      | [Optional] bowtie2 indexes of the reference genome (path with prefix)
REF_DIR           | path to your reference chromosomes directory (Ensembl format)

#### Parental files and informations:

You need to specify at least the name of the two genotypes of the parents of your hybrid strain. You also have the possibility to give your own files if you already have generated them. If optional files are left empty, files will be generated prior to mapping.

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
ID_GENO1          | name of the first genotype
FASTA_GENO1       | [Optional] genome of the first genotype (.fa)
B2_INDEX_GENO1    | [Optional] bowtie2 indexes of the first genotype
ID_GENO2          | name of the second genotype
FASTA_GENO2       | [Optional] genome of the second genotype (.fa)
B2_INDEX_GENO2    | [Optional] bowtie2 indexes of the second genotype
FASTA_DIP         | [Optional] diploid genome of the two genotypes (.fa)
B2_INDEX_DIP      | [Optional] bowtie2 indexes of the diploid genome
FASTA_NMASK       | [Optional] N-masked genome of the two genotypes (.fa)
B2_INDEX_NMASK    | [Optional] bowtie2 indexes of the N-masked genome

#### Strategy and mapper:

Here, you can select the strategy you wish to apply for the mapping of your reads and the mapper (so far, only Bowtie2 and Tophat are supported).

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
MAP_REF           | YES (1) or NO (0)
MAP_N             | 1 or 0
MAP_PAR           | 1 or 0
MAP_DIP           | 1 or 0
MAPPER            | BOWTIE2 or TOPHAT


#### Outputs:

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
GEN_OUT           | Directory where generated genome will be stored
FASTA_OUT         | sub-directory with all the generated FASTA files
INDEXES           | sub-directory with all the indexes
VCF_OUT           | sub-directory with all the generated VCF files
OUT_DIR           | [TO BE REMOVED] Output directory for alignment files

### Mapper parameters:

#### Bowtie2:

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
BOWTIE2_DIR       | path to Bowtie2 directory
B2_OPT            | Bowtie2 options

Here preset options are provided for Bowtie2 options:

```bash
B2_OPT="--reorder -p 8 --very-sensitive"
```

#### Tophat:

**VARIABLE NAME**   | **CONTENT**
------------------- | -----------
TOPHAT_DIR          | path to Tophat directory
TOPHAT_GTF          | transcript GTF file (.gtf)
TOPHAT_TRANSC_INDEX | index of transcript file
TOPHAT_OPT          | Tophat options
TOPHAT_B2_OPT       | Bowtie2 options

Here preset options are provided for Tophat and Bowtie2 options:

```bash
TOPHAT_OPT="-p 8 -g 1"
TOPHAT_B2_OPT="--b2-very-sensitive"
```

------------------------------

# How to use it ?

## Generation of genomes and indexes

Before performing mapping with any strategy, you should generate the missing files (genomes and indexes) corresponding to your `CONFIG` file with the following command:

```bash
./build_genomes_indexes.sh -c CONFIG
```

Generated files will be stored in the `GEN_OUT` directory and you will be able to run the script in parallel for this `CONFIG` file for different dataset.
*Note that you do not need to fill in the `CONFIG` file again after generation of the different files.*

## Mapping

Once the necessary files are specified or generated, mapping can be performed using the parameters below:

**VARIABLE NAME**   | **CONTENT**
------------------- | -----------
FASTQ_F             | forward reads (paired-end) or reads (single-end) (.fq/.fq.gz)
FASTQ_R             | [Paired-end only] reverse reads
OUT_DIR             | Output directory for the alignment files
OUT_NAME            | Name for the output alignment files

You simply run the following command to perform mapping:


```bash
./read_mapping.sh -f $FASTQ_F -r $FASTQ_R -o $OUT_DIR -n $OUT_NAME -c CONFIG
```

## Count table

Count table for allele-specific and total reads (common and specific) were generated with FeatureCount tool.  

```bash
./count_table.sh -b BAM_FILE -o OUT_DIR -n OUT_NAME -c CONFIG    
```

**GENE** | **PARENT1** | **PARENT2**| **ALLELIC_RATIO: PARENT1/(PARENT1+PARENT2)**|**ALL_READS**|
-------- | ---------| ---------|-------------------------|--------|
