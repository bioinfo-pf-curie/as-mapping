# Mapping Pipeline for Allele Specific analysis

This page is a quick start guide, please read the full [online manual](link) for more information.

See NEWS for information about changes in this and previous versions.

## What is it ?

This pipeline was designed for the mappings of reads prior to allele specific analysis. Most of the described methods are present in the present pipeline.

## Contact

For any question about the pipeline, please contact
<kenzo.hillion@curie.fr>

# Quick Start Guide

## Requirement

The following dependancies are required :

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 2.2.5)
* [vcf2diploid](http://alleleseq.gersteinlab.org/home.html) (version 0.2.6a)
* [samtools](http://samtools.sourceforge.net) (version 1.1)
* [clinTools](https://github.com/viv-1/clinTools)
 
## How to install it ?

To install :

```bash
tar -zxvf master.tar.gz
cd pip-master
make CONFIG_SYS=config-install.txt install
```

## CONFIG File

The `CONFIG.sh` file contains all the informations and paths for a simulation.

* The different path of the tools :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
bowtie2           | path to bowtie2 directory
vcf2diploid       | path to vcf2diploid.jar
samtools          | path to samtools
checkVariants     | path to check\_variants.py from clinTools
map\_path         | path of the current directory

* Information about the two genotypes (as the appear in the VCF file) :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
id\_geno1         | name of the first genotype
id\_geno2         | name of the second genotype

* INPUT files :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
fq\_reads         | path to your single-end reads (.fq.gz or .fq)
ref\_geno         | path to the full reference genome (.fasta)
full\_vcf         | path to the VCF file (.vcf). For mouse, VCF files are downloaded [here](ftp://ftp-mouse.sanger.ac.uk/)
ref\_dir          | path to the directory containing all reference chromosomes (.fa)


## How to use it ?

## Test dataset
