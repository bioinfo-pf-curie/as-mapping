# Single-end Reads Simulation pipeline

This page is a quick start guide, please read the full [online manual](link) for more information.

See NEWS for information about changes in this and previous versions.

## What is it ?

This pipeline was designed for the simulation of single-end reads of a heterozygous individual.

## Contact

For any question about the pipeline, please contact
<kenzo.hillion@curie.fr>

# Quick Start Guide

## How to install it ?

The following dependancies are required :

* [BEDTools](http://bedtools.readthedocs.io/en/latest/) (version 2.21.0)
* [vcf2diploid](http://alleleseq.gersteinlab.org/home.html) (version 0.2.6a)
* [samtools](http://samtools.sourceforge.net) (version 1.1)
* [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (version ART-ChocolateCherryCake-03-19-2015)
 
To install :

```bash
tar -zxvf master.tar.gz
cd pip-master
make CONFIG_SYS=config-install.txt install
```

| SYSTEM CONFIGURATION |
| -------------------- |
|PREFIX | Path to installation folder |


## Input Files

**A CONFIG file** containing:

* The different path of the tools :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
bedtools          | path to bedtools
art               | path to ART
vcf2diploid       | path to vcf2diploid.jar
samtools          | path to samtools
simreads          | path of the current directory

* Information about the two genotypes (as the appear in the VCF file) :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
id\_geno1         | name of the first genotype
id\_geno2         | name of the second genotype

* INPUT files :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
ref               | path to your reference genome (.fa)
full\_vcf         | path to the VCF file (.vcf). For mouse, VCF files are downloaded [here](ftp://ftp-mouse.sanger.ac.uk/)
mappa             | [OPTIONAL] mappability tracks (.bed). Mappability tracks can be generated using [GEMTOOLS](https://github.com/gemtools/gemtools)

* ART parameters (refer to art\_illumina\_README for more informations) :

**PARAMETER** | **DESCRIPTION**
------------- | ---------------
read\_length  | Length of the simulated reads
subs          | Correspond to the qs parameter of ART (default : 0). Influence base quality of the reads, thus substitution rate too.
rs            | Seed for random number generation
ir            | Insertion rate (default: 0.00009)
dr            | Deletion rate (default: 0.00011)i
sequencer     | Sequencing system of the built-in profile used for simulation

* OUTPUT files and directories : By default names are given but can be changed by the user

## How to use it ?

The pipeline is divided in 4 scripts allowing a modularity for the way the user wants to simulate reads.

The `biallelic_exemple.sh` script is provided as an exemple pour biallelic simulation from two strains.

## Test dataset
