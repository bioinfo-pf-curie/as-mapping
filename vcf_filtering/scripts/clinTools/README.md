# clinTools

## Requirement

This toolkit work with python 2

To use it you must install:
  * samtools version 1.1
  * hgvs python module version 0.3.7
  * pysam python module version 0.8.3

You can use the pip module manager for example in a virtual env

```
virtualenv clintools
source clintools/bin/activate
pip install hgvs==0.3.7 pysam==0.8.3
```

Don't forget to update your PATH environnement variable to use the appropriate samtools version.

## Table maker

This script and various annotation to variant files annotated with annovar.

#### Launching

This script require 5 input :

* A configuration file which is a tabulated file with 3 field per line : sample_id, annotated_file, bam_file
* Another configuration file containing an list of selected transcript NM id (one per line)
* The gtf of refGene mRNA.
* The fasta sequence of the reference genome
* A configuration file containing the association between chromosome as reference in the fasta file and the corresponding NC id.

```
cd examples/table_maker/
../../table_maker.py -Q 20 -q 20 example.conf prefered_nm.conf hg19_mRNA.gtf hg19.fa chr_accessions_hg19_GRCh37.p13.tsv > example_table_report.tsv
```


All annotated file must contain the following mandatory fields : "Start", "End", "Chr", "Otherinfo", "Ref", "Alt", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene"


#### Output file description

* Barcode : the sample_id from the firs configuration file
* Gene : Gene name from annovar "Gene.refGene" annotation
* Chr, start, end, ref , alt : informations extracted from vcf (annovar "Otherinfo" field)
* all not mandatory fields from annovar
* NM : the nm identifier of the impacted transcript
* cDNAchange : hgvs c. annotation
* AAchange : hgvs p. annotation
* all_hgvs : all hgvs annotation reported by annovar
* Status : "Heterozygous" if the allelic ratio is below 75 % and "Homozygous" otherwise.
* Depth, allelic ratio, number of reads supporting reference, alternative allele and all bases
* Strand bias for each read counts
* Grantham score : a score representing the physico-chemical effect of an amino acid substitution (Grantham 1974)
* can_splice_distance : the distance from the nearest canonic splice site
* MES_ref, MES_alt, MES_delta : maxentscan score (splicing) for reference sequence, variant sequence and ratio between 2 scores
* exon : the impacted or the nearest exon
* maxentscan donor and acceptor scores for reference sequence and variant sequence at variant position
* Base_around : the sequence around variant (useful to see stretch)
* minimum and maximum distances between variant position and start / end for reads
