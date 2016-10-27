# Single-end Reads Simulation pipeline

------------------------------

See NEWS for information about changes in this and previous versions.

## What is it ?

This pipeline was designed for the simulation of single-end reads of a heterozygous individual.

## Contact

For any question about the pipeline, please contact <kenzo.hillion@curie.fr>

------------------------------

# Quick Start Guide

## How to install it ?

#### Required dependencies:

* [BEDTools](http://bedtools.readthedocs.io/en/latest/) (version 2.21.0)
* [SNPsplit](http://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/) (version 0.3.1)
* [Samtools](http://samtools.sourceforge.net) (version 1.1)
* [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) (version ART-ChocolateCherryCake-03-19-2015)

## CONFIG File

The `CONFIG_FILE` file contains all the informations and paths for a simulation.

#### The different path of the tools :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
bedtools          | path to bedtools
art               | path to ART
samtools          | path to samtools
SNPsplit_gen      | path to SNPsplit_genome_preparation script
simreads          | path of the current directory

### Inputs:

#### Informations :

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
id_geno1          | name of the first genotype
id_geno2          | name of the second genotype
coverage          | number of reads to generate per interval

#### Files:

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
regions           | regions to generate reads from (.bed)
ref_dir           | path to your reference chromosomes directory (Ensembl format)
full_vcf          | path to the VCF file (.vcf). For mouse, VCF files are downloaded [here](ftp://ftp-mouse.sanger.ac.uk/)
mappa             | [OPTIONAL] mappability tracks (.bed). Mappability tracks can be generated using [GEMTOOLS](https://github.com/gemtools/gemtools)
ASratio           | [OPTIONAL] allelic ratio for the first genotype (.bed)


#### ART parameters (refer to `art_illumina_README` for more informations) :

**PARAMETER** | **DESCRIPTION**
------------- | ---------------
read_length   | Length of the simulated reads
subs          | Correspond to the qs parameter of ART (default : 0). Influence base quality of the reads, thus substitution rate too.
rs            | Seed for random number generation
ir            | Insertion rate (default: 0.00009)
dr            | Deletion rate (default: 0.00011)
sequencer     | Sequencing system of the built-in profile used for simulation

#### OUTPUT files and directories:

By default names are given but can be changed by the user

**PARAMETER** | **DESCRIPTION**
------------- | ---------------
main_out      | Directory where generated files (genomes, vcf...) will be stored
vcf_outdir    | sub-directory with all the generated VCF files
fasta_outdir  | sub-directory with all the generated FASTA files
tmp_outdir    | sub-directory with all the generated temporary files
art_outdir    | Directory where generated reads are stored

## Output files

Output files are generated in the `art_outdir/` directory:

* __id_geno1_id_geno2.fq.gz__: simulated single-end reads ready to be mapped.
* __id_geno1_id_geno2.bam__: alignment file like file with expected position of origin of each read.

------------------------------

# How to use it ?

First, CONFIG_FILE has to be properly filled in. Take care of using the correct genotype names as they appear in the VCF file

Once the CONFIG_FILE is correctly filled in, you run the simulation with the following command:

```bash
./main.sh -c CONFIG_FILE
```