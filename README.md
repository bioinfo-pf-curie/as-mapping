# AS-mapping

**Institut Curie - Nextflow Allele-Specific Mapping pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.8-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
<!--[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)-->

## !! UNDER CONSTRUCTION !! 

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.
The code template used for this pipeline was inspired from the [nf-core](https://nf-co.re/) project.

The goal of this pipeline is to run first level allele-specific analysis for RNA-seq or ChIP-seq data.
Two main strategies are available, either through parental mapping or through N-mask genome mapping.
This pipeline was mainly build to work on Mouse sequencing data, in conjunction with the [Mouse Genomes Project](http://www.sanger.ac.uk/science/data/mouse-genomes-project). 

### Pipline summary

1. Build allele sepcific (parental or nmask) reference genome ([`SNPsplit`](https://github.com/FelixKrueger/SNPsplit))
2. Build allele specific (parental or nmask) indexes for reads alignment ([`STAR`](https://github.com/alexdobin/STAR) / [`tophat2`](http://ccb.jhu.edu/software/tophat/index.shtml) / [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml) / ['bowtie2'](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) 
3. Align reads on reference genome ([`STAR`](https://github.com/alexdobin/STAR) / [`tophat2`](http://ccb.jhu.edu/software/tophat/index.shtml) / [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml) / ['bowtie2'](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
4. Mark duplicates ([`Picard`](https://broadinstitute.github.io/picard/))
5. Split allele specific mapped reads ([`SNPsplit`](https://github.com/FelixKrueger/SNPsplit))
6. Compute alelle specific gene counts ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
7. Generate allele specific genome track (bigWig) [`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html)

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.04.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
rnaseq v2.0.0dev
=======================================================

Usage:
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' 
nextflow run main.nf --samplePlan sample_plan --genome mm9 -profile conda


 Mandatory arguments:
   --reads [file]                Path to input data (must be surrounded with quotes)
   --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
   --genome [str]                Name of genome reference
   -profile [str]                Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

 Sequencing
   --singleEnd [bool]            Specifies that the input is single end reads

 Strandedness
   --forwardStranded [bool]      The library is forward stranded
   --reverseStranded [bool]      The library is reverse stranded
   --unStranded [bool]           The default behaviour

 Genotypes Information
   --maternal [str]              Maternal genotype (must be in the vcf file)
   --paternal [str]              Paternal genotype (must be in the vcf file)
   --nmask [bool]                Run N-mask mapping stratefy. Otherwise, parental mapping will be used
   --asfasta [file]              Skip genome preparation by specifying the allele-specific fasta file(s)
   --saveReference [bool]        Save the reference files - not done by default

 References [If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field]
   --fasta                       Path to generic reference genome 
   --vcf                         Path to vcf from Mouse Sanger Project
   --gtf                         Gene annotation (.gtf)
   --blacklist [file]            Path to black list regions (.bed).

 Mapping
   --aligner [str]               Tool for read alignments ['star', 'bowtie2', 'hisat2', 'tophat2']. Default: 'star'
   --starIndex [file]            Path to STAR index
   --bowtie2Index [file]         Path to Bowtie2 index
   --hisat2Index [file]          Path to HISAT2 index
   --tophat2Index [file]         Path to TopHat2 index

 Analysis (RNA-seq)
   --asratio                     Generate allele-specific ratio table per gene

 Analysis (ChIP-seq)
   --rmDups [bool]               Remove duplicates reads
   --bigwig                      Generate allele-specific genome-wide profile (.bigWig) 
   
   
 Other options:
   --metadata [file]             Add metadata file for multiQC report
   --outdir [file]               The output directory where the results will be saved
   -w/--work-dir [file]          The temporary directory where intermediate data will be saved
   --email [str]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
   -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

 Skip options:
   --skip_multiqc                Skip MultiQC

 =======================================================
 Available Profiles

   -profile test                Set up the test dataset
   -profile conda               Build a new conda environment before running the pipeline
   -profile toolsPath           Use the paths defined in configuration for each tool
   -profile singularity         Use the Singularity images for each process
   -profile cluster             Run the workflow on the cluster, instead of locally   
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda

```

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --genome 'mm10' --paternal '129S1_SvImJ' --maternal 'CAST_EiJ' --outdir MY_OUTPUT_DIR -profile conda

```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --genome 'mm10' --paternal '129S1_SvImJ' --maternal 'CAST_EiJ' --outdir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N asmapping"

```

### Run the pipeline step-by-step

Some of the steps as building the reference or the genome indexes can take quite a lot of time.
Therefore, the pipeline should be able to run from a reference file or pre-computed indexes.

Here are a few examples. Note that in the case of parental mapping, references or indexes can be provided using a 'comma' separatated list.

#### Skip genome(s) preparation

```
nextflow run main.nf -profile cluster,toolsPath,test --aligner 'bowtie2' --saveReference \
--asfasta '/data/tmp/asmapping_testop/CAST_EiJ_maternal_genome.fa,/data/tmp/asmapping_testop/129S1_SvImJ_paternal_genome.fa' \
--outdir /data/tmp/asmap_skipref
```

```
nextflow run main.nf -profile cluster,toolsPath,test --aligner 'bowtie2' --saveReference \
--asfasta '/data/tmp/asmapping_testop/CAST_EiJ_129S1_SvImJ_nmask_genome.fa' --nmask --outdir /data/tmp//asmap_skipref_nmask
```

#### Skip genome(s) indexing

```
nextflow run main.nf -profile cluster,toolsPath,test --aligner 'bowtie2' \
--bowtie2Index '/data/tmp/asmapping_testop/CAST_EiJ_bowtie2_index/,/data/tmp/asmapping_testop/129S1_SvImJ_bowtie2_index/' \
--outdir /data/tmp/nservant/asmap_bwt2 
```

```
nextflow run main.nf -profile cluster,toolsPath,test --aligner 'bowtie2' \
--bowtie2Index '/data/tmp/asmapping_testop/CAST_EiJ_129S1_SvImJ_bowtie2_index/' --nmask \
--outdir /data/tmp/asmap_bwt2_nmask
```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using the paths defined in the configuration for each tool (see conf.tool-path.config)
-profile toolsPath

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda

```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs.


Sample ID | Sample Name | Path R1 .fastq file | [Path R2 .fastq file]

### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (P. La Rosa, N. Servant)

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

