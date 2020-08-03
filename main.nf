#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         as-mapping
========================================================================================
 Allele-specific Mapping Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/as-mapping
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    if ("${workflow.manifest.version}" =~ /dev/ ){
       dev_mess = file("$baseDir/assets/dev_message.txt")
       log.info dev_mess.text
    }

    log.info """
    as-mapping v${workflow.manifest.version}
    ======================================================================

    Usage:
    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome mm9 -profile conda
    nextflow run main.nf --samplePlan sample_plan --genome mm9 -profile conda


    Mandatory arguments:
      --reads [file]                Path to input data (must be surrounded with quotes)
      --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
      --genome [str]                Name of genome reference
      -profile [str]                Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

    Sequencing:
      --singleEnd [bool]            Specifies that the input is single end reads

    Strandedness:
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

    Trimming
      --trimMinLen [str]            Minimum length of trimmed reads. Default: 10
      --trimQuality                 Minimum quality for reads trimming. Default: 20
 
    Mapping
      --aligner [str]               Tool for read alignments ['star', 'bowtie2', 'hisat2', 'tophat2']. Default: 'star'
      --starIndex [file]            Path to STAR index
      --bowtie2Index [file]         Path to Bowtie2 index
      --hisat2Index [file]          Path to HISAT2 index
      --tophat2Index [file]         Path to TopHat2 index

    Analysis (RNA-seq)
      --rnaseq                      Activate all RNA-seq options
      --useGtf                      Specify to use the GTF file for RNA-seq mapping
      --asratio                     Generate allele-specific ratio table per gene

    Analysis (ChIP-seq)
      --chipseq                     Activate all ChIP-seq options
      --fragmentSize                Estimated fragment length used to extend single-end reads. Default: 0
      --rmDups [bool]               Remove duplicates reads
      --bigwig                      Generate allele-specific genome-wide profile (.bigWig)
      --blacklist [file]            Path to black list regions (.bed)

    Other options:
      --metadata [file]             Add metadata file for multiQC report
      --outdir [file]               The output directory where the results will be saved
      -w/--work-dir [file]          The temporary directory where intermediate data will be saved
      --email [str]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    Skip options:
      --refOnly                     Prepare annotation only. All analysis steps are skipped
      --skipFastqc                  Skip FastQC
      --skipTrimming                Skip TrimGalore!
      --skipMultiqc                 Skip MultiQC

    =======================================================
    Available Profiles

      -profile test                Set up the test dataset
      -profile conda               Build a new conda environment before running the pipeline
      -profile toolsPath           Use the paths defined in configuration for each tool
      -profile singularity         Use the Singularity images for each process
      -profile cluster             Run the workflow on the cluster, instead of locally

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable reference genomes
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
   exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.nmask && (!params.maternal || !params.paternal)){
  exit 1, "N-mask mapping strategy requires both paternal (--paternal) and maternal (--maternal) genotypes."
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the genomes paths
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.vcf = params.genome ? params.genomes[ params.genome ].vcf ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqc_config)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * CHANNELS
 */

// Reference Genome
// Can be one or two reference genomes for nmak/parental reporting
if (params.asfasta ){
  Channel.from( params.asfasta )
         .splitCsv()
         .flatten()
         .map { file(it) }
         .into { genomeFastaStar; genomeFastaBowtie2; genomeFastaHisat2 }
}else if (params.fasta){
  Channel.fromPath("${params.fasta}")
         .ifEmpty { exit 1, "Reference Genome not found: ${params.fasta}" }
         .into { fastaParentalGenome; fastaNmaskGenome }
}


// snpFile
if ((params.starIndex || params.bowtie2Index || params.hisat2Index) && params.nmask && !params.snpFile){
  exit 1, "Using existing reference but snp file is not provided. Please define --snpFile !"       
}

if (params.snpFile){
 Channel.fromPath("${params.snpFile}")
         .ifEmpty { exit 1, "SNP file not found: ${params.snpFile}" }
         .set { chSnpFile }
}

// Genome index
// Can be one or two indexes for nmask/parental mapping
// Bowtie2/Hisat2 path are managed as string to index directory
if ( params.starIndex ){
  Channel.from( params.starIndex )
         .splitCsv()
         .flatten()
         .map { file(it) }
         .set { starIdx }
}else if ( params.bowtie2Index ){
  Channel.from( params.bowtie2Index )
         .splitCsv()
         .flatten()
         .map { file(it) }
         .set { bowtie2Idx }
}else if (params.hisat2Index ){
  Channel.from( params.hisat2Index )
         .splitCsv()
         .flatten()
         .map { file(it) }
         .set { hisat2Idx }
}

// Variant information
if (params.vcf){
  Channel.fromPath("${params.vcf}")
         .ifEmpty { exit 1, "Variant database not found: ${params.vcf}" }
         .into { vcfParentalGenome; vcfNmaskGenome }
}else{
  exit 1, "Variants (.vcf) file not found: ${params.vcf}"
}

// GTF
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtfStar; gtfHisat2Splicesites; gtfHisat2; gtfHisat2Index; gtfTophat2; chGtf }
}

// Blacklist
if (params.blacklist) {
  Channel
    .fromPath(params.blacklist, checkIfExists: true)
    .set {chBlacklistBigWig}
}else{
  chBlacklistBigWig = Channel.empty()
}

// Addition reporting information
if ( params.metadata ){
   Channel
       .fromPath( params.metadata )
       .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
       .set { ch_metadata }
}

/*
 * Create channels for input read files
 */

if (!params.refOnly){
  if(params.samplePlan){
    if(params.singleEnd){
      Channel
        .from(file("${params.samplePlan}"))
        .splitCsv(header: false)
        .map{ row -> [ row[0], [file(row[2])]] }
        .into { rawReads }
    }else{
      Channel
        .from(file("${params.samplePlan}"))
        .splitCsv(header: false)
        .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
        .into { rawReads }
    }
    params.reads=false
  }
  else if(params.readPaths){
    if(params.singleEnd){
      Channel
        .from(params.readPaths)
        .map { row -> [ row[0], [file(row[1][0])]] }
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .into { rawReads }
    } else {
      Channel
        .from(params.readPaths)
        .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .into { rawReads }
    }
  } else {
    Channel
      .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .into { rawReads }
  }

/*
 * Make sample plan if not available
 */

  if (params.samplePlan){
    chSplan = Channel.fromPath(params.samplePlan)
  }else if(params.readPaths){
    if (params.singleEnd){
      Channel
        .from(params.readPaths)
        .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
        .set{ chSplan }
    }else{
      Channel
        .from(params.readPaths)
        .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
        .set{ chSplan }
    }
  }else{
    if (params.singleEnd){
      Channel
        .fromFilePairs( params.reads, size: 1 )
        .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
         }     
         .set { chSplan }
    }else{
      Channel
        .fromFilePairs( params.reads, size: 2 )
        .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }     
        .set { chSplan }
     }
  }
}else{
  chSplan = Channel.empty()
}

// Gentoypes
if (!params.maternal && !params.paternal){
  exit 1, "Error: specify at least one genotype using --maternal or --paternal"
} 

// Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
   dev_mess = file("$baseDir/assets/dev_message.txt")
   log.info dev_mess.text
}

log.info """=======================================================

asmapping : Allele-specific Mapping workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Command Line'] = workflow.commandLine
summary['Metadata']	= params.metadata
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else if (params.reads) {
   summary['Reads']        = params.reads
}
if (! params.refOnly){
  summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
}
summary['Genome']       = params.genome
if(params.aligner == 'star'){
  summary['Aligner'] = "star"
  if(params.starIndex) summary['STAR Index'] = params.starIndex
} else if(params.aligner == 'tophat2') {
  summary['Aligner'] = "Tophat2"
  if(params.bowtie2Index) summary['Tophat2 Index'] = params.bowtie2Index
} else if(params.aligner == 'bowtie2') {
  summary['Aligner'] = "Bowtie2"
  if(params.bowtie2Index) summary['Bowtie2 Index'] = params.bowtie2Index
} else if(params.aligner == 'hisat2') {
  summary['Aligner'] = "HISAT2"
  if(params.hisat2Index) summary['HISAT2 Index'] = params.hisat2Index
}
if ( params.nmask ){
  summary['Mapping'] = 'N-mask'
}else{
  summary['Mapping'] = 'Parental'
}
summary['Paternal Genome'] = params.paternal ?: params.vcfRef
summary['Maternal Genome'] = params.maternal ?: params.vcfRef
summary['Mode'] = params.rnaseq ? "RNA-seq" : params.chipseq ? "ChIP-seq" : "None"
if ( params.reverseStranded ){
  summary['Strandness'] = 'Reverse'
}else if ( params.forwardStranded ){
  summary['Strandness'] = 'Forward'
}else{
  summary['Strandness'] = 'Unstranded'
}
summary['rmDups'] = params.rmDups || params.chipseq ? 'Yes' : 'No'
summary['asRatio'] = params.asratio || params.rnaseq ? 'Yes' : 'No'
summary['Save Reference'] = params.saveReference ? 'Yes' : params.refOnly ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * PREPROCESSING - Prepare Genome
 */

if (!params.asfasta && !params.starIndex && !params.bowtie2Index && !params.hisat2Index){

  process prepareParentalReferenceGenome {
    label 'process_highmem'
    publishDir "${params.outdir}/reference_genome", mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("_report.txt") > 0) "logs/$filename"
                 else if (params.saveReference || params.refOnly) filename
                 else null
                }
    when:
    !params.nmask

    input:
    file reference from fastaParentalGenome.collect()
    file vcf from vcfParentalGenome.collect()

    output:
    file("*_report.txt") into chGenomeParentalReport
    file("*_genome.fa") into chGenomeParentalFasta

    script:

    if (params.maternal != params.vcfRef && params.paternal != params.vcfRef){
      optsStrain = "--strain ${params.paternal} --strain2 ${params.maternal}"
    }else if (params.maternal != params.vcfRef){
      optsStrain = "--strain ${params.maternal}"
    }else if (params.paternal != params.vcfRef){
      optsStrain = "--strain ${params.paternal}"
    }
    
    if (params.paternal){
      concatPat = params.paternal == params.vcfRef ? "cat genome/*.fa > ${params.vcfRef}_paternal_genome.fa" : "cat ${params.paternal}_full_sequence/*.fa > ${params.paternal}_paternal_genome.fa"
    }
    if (params.maternal){
      concatMat = params.maternal == params.vcfRef ? "cat genome/*.fa > ${params.vcfRef}_maternal_genome.fa" : "cat ${params.maternal}_full_sequence/*.fa* > ${params.maternal}_maternal_genome.fa"
    }

    """
    mkdir genome; cd genome; faidx -x ../${reference}; cd ..
    SNPsplit_genome_preparation $optsStrain --reference_genome genome --vcf_file ${vcf} --no_nmasking
    ${concatPat}
    ${concatMat}
    """
  }

  process prepareNmaskReferenceGenome {
    label 'process_highmem'
    publishDir "${params.outdir}/reference_genome", mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("_report.txt") > 0) "logs/$filename"
                 else if (params.saveReference || params.refOnly) filename
                 else null
                }

    when:
    params.nmask

    input:
    file reference from fastaNmaskGenome.collect()
    file vcf from vcfNmaskGenome.collect()

    output:
    file("*_report.txt") into chGenomeNmaskReport
    file("*nmask.fa") into chGenomeNmaskFasta
    file("all*.txt") into chSnp

    script:
    if (params.maternal != params.vcfRef && params.paternal != params.vcfRef){
      optsStrain = "--strain ${params.paternal} --strain2 ${params.maternal}"
      nmaskPattern = "*dual_hybrid*_N-masked/*.fa"
      opref = "${params.paternal}_${params.maternal}"
     }else if (params.maternal != params.vcfRef){
      optsStrain = "--strain ${params.maternal}"
      nmaskPattern = "*_N-masked/*.fa"
      opref = "${params.vcfRef}_${params.maternal}"
     }else if (params.paternal != params.vcfRef){
      optsStrain = "--strain ${params.paternal}"
      nmaskPattern = "*_N-masked/*.fa"
      opref = "${params.paternal}_${params.vcfRef}"
     }

    """
    mkdir genome; cd genome; faidx -x ../${reference}; cd ..
    SNPsplit_genome_preparation $optsStrain --reference_genome genome --vcf_file ${vcf} -nmasking
    cat ${nmaskPattern} > ${opref}_nmask.fa
    gunzip all*gz
    """                                            
  }

  if (params.nmask){
    chGenomeNmaskFasta.into{genomeFastaStar; genomeFastaBowtie2; genomeFastaHisat2}
    chGenomeReport = chGenomeNmaskReport
    chSnpFile = chSnp
  }else{
    chGenomeParentalFasta.into{genomeFastaStar; genomeFastaBowtie2; genomeFastaHisat2}
    chGenomeReport = chGenomeParentalReport
    chSnpFile = Channel.empty()
  }
}else{
  if(!params.nmask){
    chSnpFile = Channel.empty()
  }
}

/*
 * PREPROCESSING - Build Index
 */ 

if ( params.aligner == 'star' && !params.starIndex ){
  process makeStarIndex {
    label 'process_high'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
       saveAs: {filename -> if (params.saveReference || params.refOnly) filename else null }

    input:
    file(fasta) from genomeFastaStar.flatten()

    output:
    file "*STAR_index" into starIdx

    script:
    strainPrefix = fasta.toString() - ~/(\_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
    """
    mkdir -p ${strainPrefix}_STAR_index
    STAR --runMode genomeGenerate --limitGenomeGenerateRAM 33524399488 --runThreadN ${task.cpus} --genomeDir ${strainPrefix}_STAR_index --genomeFastaFiles $fasta
    """
  }
}

if ( (params.aligner == "bowtie2" || params.aligner == "tophat2") && !params.bowtie2Index ){
  process makeBowtie2Index {
    label 'process_medium'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
       saveAs: {filename -> if (params.saveReference || params.refOnly) filename else null }

    input:
    file(fasta) from genomeFastaBowtie2.flatten()

    output:
    file("${strainPrefix}_bowtie2_index") into (bowtie2Idx, tophat2Idx)

    script:
    strainPrefix = fasta.toString() - ~/(\_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
    base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir -p ${strainPrefix}_bowtie2_index
    bowtie2-build ${fasta} ${strainPrefix}_bowtie2_index/${base}
    """
  }
}

if ( params.aligner == 'hisat2' && !params.hisat2Index ){
  process makeHisat2Splicesites {
    label 'process_low'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
       saveAs: {filename -> if (params.saveReference || params.refOnly) filename else null }

    input:
    file gtf from gtfHisat2Splicesites.collect()

    output:
    file "${gtf.baseName}.hisat2_splice_sites.txt" into indexingSplicesites, alignmentSplicesites

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
    """
  }

  process makeHisat2Index {
    label 'process_extra'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
         saveAs: {filename -> if (params.saveReference || params.refOnly) filename else null }

    input:
    file fasta from genomeFastaHisat2.flatten()
    file indexing_splicesites from indexingSplicesites.collect()
    file gtf from gtfHisat2Index.collect()

    output:
    file("${strainPrefix}_hisat2_index") into hisat2Idx

    script:
    strainPrefix = fasta.toString() - ~/(\_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
    base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir -p ${strainPrefix}_hisat2_index
    hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt
    hisat2-build -p ${task.cpus} --ss $indexing_splicesites --exon ${gtf.baseName}.hisat2_exons.txt $fasta ${strainPrefix}_hisat2_index/${base}
    """
  }
}

/*
 * TRIMMING
 */

process trimGalore {
  tag "$name" 
  publishDir "${params.outdir}/trimming", mode: 'copy',
              saveAs: {filename -> filename.indexOf(".log") > 0 ? "logs/$filename" : "$filename"}
  when:
  !params.skipTrimming

  input:
  set val(name), file(reads) from rawReads

  output:
  set val(name), file("*fastq.gz") into readsFastqc, readsStar, readsBowtie2, readsHisat2, readsTophat2
  set val(name), file("*trimming_report.txt") into trimResults

  script:
  prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(\.fq)?(\.fastq)?(\.gz)?$/
  if (params.singleEnd) {
    """
    trim_galore --trim-n \
                --quality ${params.trimQuality} \
                --length ${params.trimMinlen} \
                --gzip $reads --basename ${prefix} --cores ${task.cpus}
    mv ${prefix}_trimmed.fq.gz ${prefix}_trimmed_R1.fastq.gz
    """

  }else {
    """
    trim_galore --trim-n \
                --quality ${params.trimQuality} \
                --length ${params.trimMinLen} \
                --paired --gzip $reads --basename ${prefix} --cores ${task.cpus}
    mv ${prefix}_R1_val_1.fq.gz ${prefix}_trimmed_R1.fastq.gz
    mv ${prefix}_R2_val_2.fq.gz ${prefix}_trimmed_R2.fastq.gz
    """
  }
}

if (params.skipTrimming){
  rawReads.into{ readsFastqc, readsStar, readsBowtie2, readsHisat2, readsTophat2 }
}


/*
 * FastQC
 */

process fastqc {
    tag "$name (raw)"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skipFastqc

    input:
    set val(name), file(reads) from readsFastqc

    output:
    file "*_fastqc.{zip,html}" into fastqcResults

    script:
    """
    fastqc -q $reads -t ${task.cpus}
    """
}


/*
 * MAPPING
 */

// STAR
if ( params.aligner == 'star' && !params.refOnly ){
  process star {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else filename}

    input:
    set val(prefix), file(reads) from readsStar
    each index from starIdx
    file gtf from gtfStar.collect().ifEmpty([])

    output:
    set val(prefix), file ("*bam") into starBam
    file "*Log.final.out" into starLogs

    script:
    def gtfOpts = params.gtf && (params.useGtf || params.rnaseq) ? "--sjdbGTFfile $gtf" : ''
    def mandatoryOpts = "--alignEndsType EndToEnd --outSAMattributes NH HI NM AS MD --outSAMtype BAM Unsorted --outSAMunmapped Within"
    def genomeBase = index.baseName - ~/_STAR_index$/
    """
    STAR --genomeDir $index \\
       ${gtfOpts} \\
       ${mandatoryOpts} \\
       --readFilesIn $reads  \\
       --runThreadN ${task.cpus} \\
       --runMode alignReads \\
       --readFilesCommand zcat \\
       --runDirPerm All_RWX \\
       --outTmpDir /local/scratch/asmap_\$(date +%d%s%S%N) \\
       --outFileNamePrefix ${prefix}_${genomeBase}  \\
       --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina

    ## sort by read name
    samtools sort -n ${prefix}_${genomeBase}Aligned.out.bam -@ ${task.cpus} -T ${prefix} -o ${prefix}_${genomeBase}_sorted.bam 
    mv ${prefix}_${genomeBase}_sorted.bam ${prefix}_${genomeBase}Aligned.out.bam
    """
  }
}

// Bowtie2
if ( params.aligner == 'bowtie2'  && !params.refOnly){
  process bowtie2 {
    tag "$prefix"
    label 'process_medium'
    publishDir "${params.outdir}/alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else filename}

    input:
    set val(prefix), file(reads) from readsBowtie2
    each file(index) from bowtie2Idx

    output:
    set val(prefix), file("*.bam") into bowtie2Bam
    file("*.log") into bowtie2Logs

    script:
    idxDir = file(index.toRealPath())
    allFiles = idxDir.list()
    genomeBase = allFiles[0] - ~/(\.rev)?(.\d.bt2)/ 
    bwt2Opts = params.nmask ? "-D 70 -R 3 -N 0 -L 20 -i S,1,0.50" : ""
    inCmd = params.singleEnd ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 --very-sensitive --end-to-end --reorder \\
            ${bwt2Opts} \\
            --rg-id BMG --rg SM:${prefix} \\
            -p ${task.cpus} \\
            -x ${index}/${genomeBase} \\
            ${inCmd} 2> ${prefix}_${genomeBase}_bowtie2.log | samtools view -bS - > ${prefix}_${genomeBase}.bam
    """
  }
}


// TopHat2
if(params.aligner == 'tophat2' && !params.refOnly){
  process tophat2 {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".align_summary.txt") > 0) "logs/$filename"
            else filename
        }

    input:
    set val(prefix), file(reads) from readsTophat2 
    each file(index) from tophat2Idx
    file gtf from gtfTophat2.collect()

    output:
    set val(prefix), file("${prefix}_${genomeBase}.bam") into tophat2Bam
    file "*.align_summary.txt" into tophat2Logs

    script:
    idxDir = file(index.toRealPath())
    allFiles = idxDir.list()
    genomeBase = allFiles[0] - ~/(\.rev)?(.\d.bt2)/ 

    def stranded_opt = '--library-type fr-unstranded'
    if (params.forwardStranded){
      stranded_opt = '--library-type fr-secondstrand'
    }else if (params.reverseStranded){
      stranded_opt = '--library-type fr-firststrand'
    }
    def out = './mapping'
    def sample = "--rg-id ${prefix} --rg-sample ${prefix} --rg-library Illumina --rg-platform Illumina --rg-platform-unit ${prefix}"
    """
    mkdir -p ${out}
    export PATH=/bioinfo/local/build/Centos/python/python-2.7.13/bin/:/bioinfo/local/build/Centos/tophat/tophat_2.1.1/bin/:/bioinfo/local/build/Centos/bowtie2/bowtie2-2.2.9/:$PATH
    tophat2 -p ${task.cpus} \\
             ${sample} \\
             --GTF $gtf \\
             ${stranded_opt} \\
             -o ${out} \\
             ${index}/${genomeBase} \\
             ${reads} 

    mv ${out}/accepted_hits.bam ${prefix}_${genomeBase}.bam
    mv ${out}/align_summary.txt ${prefix}_${genomeBase}.align_summary.txt
    """
  }
}


// HiSat2
if ( params.aligner == 'hisat2' && !params.refOnly){
  process hisat2 {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
	        else filename}
    input:
    set val(prefix), file(reads) from readsHisat2
    each file(index) from hisat2Idx
    file alignmentSplicesites from alignmentSplicesites.collect()

    output:
    set val(prefix), file("${prefix}_${genomeBase}.bam") into hisat2Bam
    file "*hisat2_summary.txt" into hisat2Logs

    script:
    idxDir = file(index.toRealPath())
    allFiles = idxDir.list()
    genomeBase = allFiles[0] - ~/(\.rev)?(.\d.ht2)/
    inCmd = params.singleEnd ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def strandness = ''
    if (params.forwardStranded){
        strandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (params.reverseStranded){
        strandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    """
    hisat2 -x $index/${genomeBase} \\
           $inCmd \\
           $strandness \\
           --known-splicesite-infile $alignmentSplicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
	   --no-softclip \\
           --summary-file ${prefix}_${genomeBase}.hisat2_summary.txt \\
           | samtools view -bS -F 4 -F 256 - > ${prefix}_${genomeBase}.bam
    """
  }
}

if (params.refOnly){
  chNmaskBams = Channel.empty()
  chBams = Channel.empty()
  chMappingMqc = Channel.empty()
}else{
  if (params.aligner == 'star'){
    starBam.into{chNmaskBams; chBams}
    chMappingMqc = starLogs
  }else if (params.aligner == 'bowtie2'){
    bowtie2Bam.into{chNmaskBams; chBams}
    chMappingMqc = bowtie2Logs
  }else if (params.aligner == 'hisat2'){
    hisat2Bam.into{chNmaskBams; chBams}
    chMappingMqc = hisat2Logs
  }else if (params.aligner == 'tophat2'){
    tophat2Bam.into{chNmaskBams; chBams}
    chMappingMqc = tophat2Logs
  }
}

/*****************************
 * Tag allele-specific reads
 */

// Parental mapping
process tagParentalBams {
  tag "${prefix}"
  label 'process_low'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename}

  when:
  !params.nmask

  input:
  set val(prefix), file(bams) from chBams.groupTuple()

  output:
  set val(prefix), file("${prefix}_parentalMerged.bam") into chTagParentalBams
  file("*.log") into tagParentalLog

  script:
  pbam = bams[0] =~ params.paternal ? bams[0] : bams[1]
  mbam = bams[0] =~ params.maternal ? bams[0] : bams[1]
  """
  mergeParentalMapping.py -m ${mbam} -p ${pbam} -o ${prefix}_parentalMerged.bam
  mv mergeAlignReport.log ${prefix}_mergeAlignReport.log
  """
}


// N-mask mapping
process tagNmaskBams {
  tag "${prefix}"
  label 'process_mediummem'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
               saveAs: {filename ->
               if (filename.indexOf(".txt") > 0) "logs/$filename"
               else filename}

  when:
  params.nmask

  input:
  set val(prefix), file(bam) from chNmaskBams
  file(snpFile) from chSnpFile.collect()

  output:
  set val(prefix), file("*allele_flagged.bam") into chTagNmaskBams
  file("*report.txt") into tagNmaskLog

  script:
  opts = params.singleEnd ? "" : "--paired --singletons"
  """
  SNPsplit ${opts} --snp_file ${snpFile} ${bam}
  """
}

if (params.nmask){
  chTagNmaskBams.into{chTagBams ; chTagBamsPicard}
  tagNmaskLog.set{chTagLog}
}else{
  chTagParentalBams.into{chTagBams ; chTagBamsPicard}
  tagParentalLog.set{chTagLog}
}


/****************************
 * BAM filering
 */

process rmDuplicates{
  tag "${prefix}"
  label 'process_medium'
  publishDir path: "${params.outdir}/alignment", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && !filename.endsWith(".bam.bai") ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) ) filename
             else null}

  when:
  params.rmDups || params.chipseq

  input:
  set val(prefix), file(bams) from chTagBamsPicard

  output:
  set val(prefix), file("*marked.bam") into chMarkedBams
  file "*.txt" into chMarkedPicstats

  script:
  pfix=bams.toString() - ~/(.bam)?$/
  """
  ## Sort by coordinates
  samtools sort ${bams} -@ ${task.cpus} -T ${prefix} -o ${pfix}_sorted.bam
  picard -Xmx4g MarkDuplicates \\
    INPUT=${pfix}_sorted.bam \\
    OUTPUT=${pfix}_marked.bam \\
    REMOVE_DUPLICATES=true \\
    METRICS_FILE=${pfix}.MarkDuplicates.metrics.txt \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR=tmpdir

  ## Back to query name sort
  samtools sort -n ${pfix}_marked.bam -@ ${task.cpus} -T ${prefix} -o ${pfix}_marked_sorted.bam
  mv ${pfix}_marked_sorted.bam ${pfix}_marked.bam
  """
}

if (params.rmDups || params.chipseq){
  chMarkedBams.into{chFiltBams; chFiltBamsSplit}
}else{
  chTagBams.into{chFiltBams; chFiltBamsSplit}
  chMarkedPicstats = Channel.empty()
}


/***************************
 * Split genome1/genome2
 */ 
process splitTaggedBam {
  tag "${prefix}"
  label 'process_low'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename}

  input:
  set val(prefix), file(asBam) from chFiltBamsSplit

  output:
  set val(prefix), file("*genome1.bam"), file("*genome2.bam") into genomeBams
  file("*sort.txt") into splitLogs 
  file("*sort.mqc") into splitMqc 

  script:
  opts = params.singleEnd ? "" : "--paired --singletons"
  """
  tag2sort ${opts} ${asBam}
  snpsplit2mqc.sh *sort.txt > ${prefix}_sort.mqc

  """
}

genomeBams.join(chFiltBams).into{chBamCount; chBamWig}


/*************************
 * AS ratio
 */

process geneASratio {
  tag "${prefix}"
  label 'process_medium'  
  publishDir "${params.outdir}/asratio", mode: 'copy',

  when:
  params.asratio || params.rnaseq

  input:
  set val(prefix), file(bam1), file(bam2), file(bamTag) from chBamCount
  file(gtf) from chGtf.collect()

  output:
  file("*count.txt") into asCounts
  file("*allelicRatio.txt") into asRatio
  file("*average.mqc") into asRatioMqc

  script:
  opts=params.singleEnd ? "" : "-p"
  strandness = "-s 0"
  if (params.forwardStranded){
      strandness = "-s 1"
  } else if (params.reverseStranded){
      strandness = "-s 2"
  }
  """
  featureCounts ${opts} ${strandness} -T ${task.cpus} -a ${gtf} -g 'gene_name' -o ${prefix}_count.txt ${bam1} ${bam2} ${bamTag}
  awk -F '\\t' 'BEGIN{OFS="\t"; print "Gene", "Genome1", "Genome2", "ASratio", "allReads"}{\
      if(\$1!~/^#/ && \$1!="Geneid"){if(\$7+\$8>0){as=\$7/(\$7+\$8)} else{as="NA"};\
      print \$1,\$7,\$8,as,\$9}}' ${prefix}_count.txt >  ${prefix}_allelicRatio.txt
  awk -F'\\t' '(NR>2 && (\$7+\$8)>0){n=split(\$2,chrom,";");as=\$7/(\$7+\$8);mr[chrom[1]]+=as;cr[chrom[1]]+=1}\
      END{for(x in mr){print x","mr[x]/cr[x]}}' ${prefix}_count.txt | sort -k1,1V >  ${prefix}_average.mqc
  """
}


/***********************
 * UCSC tracks
 */

process bigWig {
  tag "${prefix}"
  label 'process_medium'
  publishDir "${params.outdir}/bigWig", mode: "copy"

  when:
  params.bigwig || params.chipseq

  input:
  set val(prefix), file(bam1), file(bam2), file(bamTag) from chBamWig
  file(BLbed) from chBlacklistBigWig.collect().ifEmpty([])

  output:
  set val(prefix), file('*.bigwig') into chBigWig

  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = "--extendReads"
  }    
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  """
  nbreads=\$(samtools view -c ${bamTag})
  sf=\$(echo "10000000 \$nbreads" | awk '{printf "%.2f", \$1/\$2}')
  
  samtools sort -@ ${task.cpus} -T ${prefix} -o ${prefix}_genome1_sorted.bam ${bam1}
  samtools index ${prefix}_genome1_sorted.bam
  bamCoverage -b ${prefix}_genome1_sorted.bam \\
              -o ${prefix}_genome1_norm.bigwig \\
              -p ${task.cpus} \\
	      ${extend} \\
              ${blacklistParams} \\
              --scaleFactor \$sf

  samtools sort -@ ${task.cpus} -T ${prefix} -o ${prefix}_genome2_sorted.bam ${bam2}
  samtools index ${prefix}_genome2_sorted.bam
  bamCoverage -b ${prefix}_genome2_sorted.bam \\
              -o ${prefix}_genome2_norm.bigwig \\
              -p ${task.cpus} \\
	      ${extend} \\
              ${blacklistParams} \\
              --scaleFactor \$sf
  """
}


/*****************************************************
 * MultiQC
 */


process get_software_versions {
  output:
  file 'software_versions_mqc.yaml' into software_versions_yaml

  script:
  """
  echo $workflow.manifest.version &> v_main.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  hisat2 --version &> v_hisat2.txt
  #tophat2 --version &> v_tophat2.txt
  echo "TopHat vXXX" > v_tophat2.txt
  STAR --version &> v_star.txt
  bowtie2 --version &> v_bowtie2.txt
  picard MarkDuplicates --version &> v_markduplicates.txt || true
  echo \$(bamCoverage --version 2>&1) > v_deeptools.txt || true
  featureCounts -v &> v_featurecounts.txt
  samtools --version &> v_samtools.txt
  multiqc --version &> v_multiqc.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflow_summary_mqc {
  when:
  !params.skipMultiqc

  output:
  file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/as-mapping'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
  """.stripIndent()
}

process multiqc {
    label 'process_low'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skipMultiqc

    input:
    file splan from chSplan.collect()
    file metadata from ch_metadata.ifEmpty([])
    file multiqc_config from chMultiqcConfig    
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()
    file ('trimming/*') from trimResults.collect().ifEmpty([])
    file ('alignment/*') from chMappingMqc.collect().ifEmpty([])
    file ('alignment/*') from chMarkedPicstats.collect().ifEmpty([])
    file ('tag/*') from chTagLog.collect().ifEmpty([])
    file ('snpsplit/*') from splitMqc.collect().ifEmpty([])    
    file ('asratio/*') from asRatioMqc.collect().ifEmpty([])

    output:
    file splan
    file "*report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName + "_asmap_report" : "--filename asmap_report"
    metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
    nmaskOpt = params.nmask ? "-n" : ""
    peOpt  = params.singleEnd ? "" : "-p"
    strandness = params.forwardStranded ? 'forward' : params.reverseStranded ? 'reverse' : 'unstranded'
    """
    stats2multiqc.sh -s ${splan} -a ${params.aligner} -d ${strandness} -f ${params.paternal} -m ${params.maternal} ${nmaskOpt} ${peOpt}
    mqc_header.py --splan ${splan} --name "Allele-Specific Mapping" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
    multiqc . -f $rtitle $rfilename -c $multiqc_config -c multiqc-config-header.yaml
    """    
}


/*
 * Sub-routine
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from chOutputDocs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

workflow.onComplete {

    /*pipeline_report.html*/

    def report_fields = [:]
    report_fields['version'] = workflow.manifest.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/oncomplete_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()
    
    // Render the HTML template
    def hf = new File("$baseDir/assets/oncomplete_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << report_txt }

    /*oncomplete file*/

    File woc = new File("${params.outdir}/workflow.oncomplete.txt")
    Map endSummary = [:]
    endSummary['Completed on'] = workflow.complete
    endSummary['Duration']     = workflow.duration
    endSummary['Success']      = workflow.success
    endSummary['exit status']  = workflow.exitStatus
    endSummary['Error report'] = workflow.errorReport ?: '-'
    String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
    println endWfSummary
    String execInfo = "${fullSum}\nExecution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
    woc.write(execInfo)

    /*final logs*/
    if(workflow.success){
        log.info "[ASmapping] Pipeline Complete"
    }else{
        log.info "[ASmapping] FAILED: $workflow.runName"
        if( workflow.profile == 'test'){
            log.error "====================================================\n" +
                    "  WARNING! You are running with the profile 'test' only\n" +
                    "  pipeline config profile, which runs on the head node\n" +
                    "  and assumes all software is on the PATH.\n" +
                    "  This is probably why everything broke.\n" +
                    "  Please use `-profile test,conda` or `-profile test,singularity` to run on local.\n" +
                    "  Please use `-profile test,conda,cluster` or `-profile test,singularity,cluster` to run on your cluster.\n" +
                    "============================================================"
        }
    }
 
}
