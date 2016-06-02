# Quick Start Guide for single-end reads Simulation

This page is a quick start guide, please read the full [online manual](link) for more information.

See NEWS for information about changes in this and previous versions.

## What is it ?

This pipeline was designed for the simulation of single-end reads.

For any question about the pipeline, please contact
<kenzo.hillion@curie.fr>

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

* The different path of the tools

```bash
bedtools=path_to_bedtools
art=path_to_desired_art
vcf2diploid=path_to_vcf2diploid.jar
samtools=path_to_samtools
extract_SNPs=path_to_extract_snps_forNmask.py
simreads=directory_path
```

## How to use it ?

First have a look at the help message !

```
command --help
usage : command -h
Use option -h|--help for more information
```

## Test dataset
