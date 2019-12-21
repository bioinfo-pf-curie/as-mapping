# AS-mapping

**Institut Curie - Nextflow Allele-Specific Mapping pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.7-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
<!--[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)-->

## !! UNDER CONSTRUCTION !! 

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

### Pipline summary

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.04.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
rnaseq v2.0.0dev
=======================================================

Usage:
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' 
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
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --genome 'hg19' --outdir MY_OUTPUT_DIR -profile conda

```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --genome 'hg19' --outdir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N rnaseq-2.0

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

