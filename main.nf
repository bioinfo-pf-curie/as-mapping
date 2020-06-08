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

    References:
      --maternal [str]
      --paternal [str]
      --nmask [bool]
      --asfasta [file]
      --saveReference [bool]        Save the reference files - not done by default

    Mapping:
      --aligner [str]               Tool for read alignments ['star', 'bowtie2', 'hisat2', 'tophat2']. Default: 'star'
      --starIndex [file]            Path to STAR index
      --bowtie2Index [file]         Path to Bowtie2 index
      --hisat2Index [file]          Path to HISAT2 index
      --tophat2Index [file]         Path to TopHat2 index

    Analysis
      --asratio                     Generate allele-specific ratio table per gene
      --bigwig                      Generate allele-specific genome-wide profile (.bigWig)
      --blacklist [file]            Path to black list regions (.bed).


    Other options:
      --rmDups                      TODO
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
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

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
         .into { fastaGenomeParental; fastaGenomeNmask }
}

if (params.snpFile){
 Channel.fromPath("${params.snpFile}")
         .ifEmpty { exit 1, "SNP file not found: ${params.snpFile}" }
         .into { chSnpFile }
}else{
  chSnpFile=Channel.empty()
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
         .into { vcfGenomeParental; vcfGenomeNmask }
}else{
  exit 1, "Variants (.vcf) file not found: ${params.vcf}"
}

// GTF
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtfStarIndex; gtfStar; gtfHisat2Splicesites; gtfHisat2; gtfHisat2Index; gtfTophat2; chGtf }
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

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         .into { rawReadsStar; rawReadsHisat2; rawReadsBowtie2; rawReadsTophat2 }
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .into { rawReadsStar; rawReadsHisat2; rawReadsBowtie2; rawReadsTophat2 }
   }
   params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { rawReadsStar; rawReadsHisat2; rawReadsBowtie2; rawReadsTophat2 }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { rawReadsStar; rawReadsHisat2; rawReadsBowtie2; rawReadsTophat2 }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { rawReadsStar; rawReadsHisat2; rawReadsBowtie2; rawReadsTophat2 }
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  ch_splan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .set{ ch_splan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ ch_splan }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { ch_splan }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { ch_splan }
   }
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
}else{
   summary['Reads']        = params.reads
}
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
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
  summary['Mapping Strategy'] = 'N-mask'
}else{
  summary['Mapping Strategy'] = 'Parental'
}
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
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
  process prepareReferenceGenome {
    label 'process_highmem'
    publishDir "${params.outdir}/reference_genome", mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("_report.txt") > 0) "logs/$filename"
                 else if (params.saveReference) filename
                 else null
                }

    input:
    file reference from fastaGenomeParental.collect()
    file vcf from vcfGenomeParental.collect()

    output:
    file("*_report.txt") into parentalGenomeReport
    file("*genome.fa") into (genomeFastaStar, genomeFastaBowtie2, genomeFastaHisat2)
    file("all*.txt") into chSnpFile  

    script:
    if (params.maternal && params.paternal){
      opts_strain = "--strain ${params.paternal} --strain2 ${params.maternal}"
    }else{
      opts_strain = params.maternal ? " --strain ${params.maternal}" : "--strain ${params.paternal}" 
    }
    if (!params.nmask){
    """
    SNPsplit_genome_preparation $opts_strain --reference_genome ${reference} --vcf_file ${vcf} --no_nmasking
    cat ${params.paternal}_full_sequence/*.fa > ${params.paternal}_paternal_genome.fa
    cat ${params.maternal}_full_sequence/*.fa > ${params.maternal}_maternal_genome.fa
    """
    }else{
    """
    SNPsplit_genome_preparation $opts_strain --reference_genome ${reference} --vcf_file ${vcf} -nmasking
    cat *dual_hybrid*_N-masked/*.fa > ${params.maternal}_${params.paternal}_nmask_genome.fa
    """
    }
  }
}

/*
 * PREPROCESSING - Build Index
 */ 

if ( params.aligner == 'star' && !params.starIndex ){
  process makeStarIndex {
    label 'process_high'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
       saveAs: {filename -> if (params.saveReference) filename else null }

    input:
    file(fasta) from genomeFastaStar.flatten()
    file gtf from gtfStarIndex.collect()

    output:
    file "*STAR_index" into starIdx

    script:
    strainPrefix = fasta.toString() - ~/(\_nmask_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
    """
    mkdir -p ${strainPrefix}_STAR_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir ${strainPrefix}_STAR_index --genomeFastaFiles $fasta
    """
  }
}

if ( (params.aligner == "bowtie2" || params.aligner == "tophat2") && !params.bowtie2Index ){
  process makeBowtie2Index {
    label 'process_low'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
       saveAs: {filename -> if (params.saveReference) filename else null }

    input:
    file(fasta) from genomeFastaBowtie2.flatten()

    output:
    file("${strainPrefix}_bowtie2_index") into (bowtie2Idx, tophat2Idx)

    script:
    strainPrefix = fasta.toString() - ~/(\_nmask_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
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
       saveAs: {filename -> if (params.saveReference) filename else null }

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
    label 'process_high'
    publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
         saveAs: {filename -> if (params.saveReference) filename else null }

    input:
    file fasta from genomeFastaHisat2.flatten()
    file indexing_splicesites from indexingSplicesites.collect()
    file gtf from gtfHisat2Index.collect()

    output:
    file("${strainPrefix}_hisat2_index") into hisat2Idx

    script:
    strainPrefix = fasta.toString() - ~/(\_nmask_genome.fa)?(\_paternal_genome.fa)?(\_maternal_genome.fa)?$/
    base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir -p ${strainPrefix}_hisat2_index
    hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt
    hisat2-build -p ${task.cpus} --ss $indexing_splicesites --exon ${gtf.baseName}.hisat2_exons.txt $fasta ${strainPrefix}_hisat2_index/${base}
    """
  }
}


/*
 * MAPPING
 */

// STAR
if ( params.aligner == 'star' ){
  process star {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/mapping", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else filename}

    input:
    set val(prefix), file(reads) from rawReadsStar
    each index from starIdx
    file gtf from gtfStar.collect().ifEmpty([])

    output:
    set val(prefix), file ('*.bam') into starBam
    file "*Log.final.out" into starLogs

    script:
    def gtfOpts = params.gtf ? "--sjdbGTFfile $gtf" : ''
    def mandatoryOpts = "--alignEndsType EndToEnd --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted"
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
       --outTmpDir /local/scratch/rnaseq_\$(date +%d%s%S%N) \\
       --outFileNamePrefix ${prefix}_${genomeBase}  \\
       --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
    """
  }
}

// Bowtie2
if ( params.aligner == 'bowtie2' ){
  process bowtie2 {
    tag "$prefix"
    label 'process_medium'
    publishDir "${params.outdir}/mapping", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else filename}

    input:
    set val(prefix), file(reads) from rawReadsBowtie2
    each file(index) from bowtie2Idx

    output:
    set val(prefix), file("*.bam") into bowtie2Bam

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
            ${inCmd} | samtools view -bS - > ${prefix}_${genomeBase}.bam
    """
  }
}


// TopHat2
if(params.aligner == 'tophat2'){
  process tophat2 {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/mapping", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".align_summary.txt") > 0) "logs/$filename"
            else filename
        }

    input:
    set val(prefix), file(reads) from rawReadsTophat2 
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
if ( params.aligner == 'hisat2' ){
  process hisat2 {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/mapping", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
	        else filename}
    input:
    set val(prefix), file(reads) from rawReadsHisat2
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

if (params.nmask){
  chBams = Channel.empty()
  if (params.aligner == 'star'){
    chNmaskBams = starBam
  }else if (params.aligner == 'bowtie2'){
    chNmaskBams = bowtie2Bam
  }else if (params.aligner == 'hisat2'){
    chNmaskBams = hisat2Bam
  }else if (params.aligner == 'tophat2'){
    chNmaskBams = tophat2Bam
 }
}else{
  chNmaskBams = Channel.empty()
  if (params.aligner == 'star'){
    chBams = starBam.groupTuple()
  }else if (params.aligner == 'bowtie2'){
    chBams = bowtie2Bam.groupTuple()
  }else if (params.aligner == 'hisat2'){
    chBams = hisat2Bam.groupTuple()
  }else if (params.aligner == 'tophat2'){
    chBams = tophat2Bam.groupTuple()
  }
}

/*****************************
 * Tag allele-specific reads
 */

process mergeParentalBams {
  tag "${prefix}"
  label 'process_low'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename}

  when:
  !params.nmask

  input:
  set val(prefix), file(bams) from chBams  

  output:
  set val(prefix), file("${prefix}_parentalMerged.bam") into parentalBams, parentalBams2split
  file("*.log") into parentalLog

  script:
  opts = params.singleEnd ? "" : "--pairedEnd"
  """
  mergeParentalMapping.py -p ${bams[0]} -m ${bams[1]} -o ${prefix}_parentalMerged.bam $opts
  mv mergeAlignReport.log ${prefix}_mergeAlignReport.log
  """
}
                                                                                                                                                                                                     
process splitTaggedBam {
  tag "${prefix}"
  label 'process_low'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename}

  when:
  !params.nmask

  input:
  set val(prefix), file(asBam) from parentalBams2split

  output:
  set val(prefix), file("*genome1.bam"), file("*genome2.bam") into genomeParentalBams

  script:
  opts = params.singleEnd ? "" : "--paired --singletons"
  """
  tag2sort ${opts} ${asBam}
  """
}                 

process snpSplitTag {
  tag "${prefix}"
  label 'process_mediummem'
  publishDir "${params.outdir}/taggedBam", mode: 'copy',
               saveAs: {filename ->
               if (filename.indexOf(".log") > 0) "logs/$filename"
               else filename}

  when:
  params.nmask

  input:
  set val(prefix), file(bam) from chNmaskBams
  file(snpFile) from chSnpFile.collect()

  output:
  set val(prefix), file("*allele_flagged.bam") into nmaskBams
  set val(prefix), file("*genome1.bam"), file("*genome2.bam") into genomeNmaskBams

  script:
  opts = params.singleEnd ? "" : "--paired --singletons"
  """
  SNPsplit ${opts} --snp_file ${snpFile} ${bam}
  """
}

if (params.nmask){
  genomeNmaskBams.join(nmaskBams).into{chBamCount; chBamWig}
}else{
  genomeParentalBams.join(parentalBams).into{chBamCount; chBamWig}
}


/*************************
 * AS ratio
 */

process geneASratio {
  tag "${prefix}"
  label 'process_medium'  

  when:
  params.asratio

  input:
  set val(prefix), file(bam1), file(bam2), file(bamTag) from chBamCount
  file(gtf) from chGtf.collect()

  output:
  file("*count.txt") into asCounts
  file("*allelicRatio.txt") into asRatio

  script:
  opts=params.singleEnd ? "" : "-C -p -P"
  strandness = "-s 0"
  if (params.forwardStranded){
      strandness = "-s 1"
  } else if (params.reverseStranded){
      strandness = "-s 2"
  }
  """
  featureCounts ${opts} ${strandness} -T ${task.cpus} -a ${gtf} -o ${prefix}_count.txt ${bam1} ${bam2} ${bamTag}
  awk -F '\\t' -v threshold=2 'BEGIN{OFS="\t"; print "Gene", "Genome1", "Genome2", "allelicRatio_Genome1/all", "allReads"}{\
      if(\$1!~/^#/ && \$1!="Geneid"){if(\$7+\$8>0 && \$7+\$8>=threshold){as=\$7/(\$7+\$8)}\
      else{as="NA"};print \$1,\$7,\$8,as,\$9}}' ${prefix}_count.txt >  ${prefix}_allelicRatio.txt
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
  params.bigwig

  input:
  set val(prefix), file(bam1), file(bam2), file(bamTag) from chBamWig
  file(BLbed) from chBlacklistBigWig.collect().ifEmpty([])

  output:
  set val(prefix), file('*.bigwig') into chBigWig

  script:
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  """
  nbreads=\$(samtools view -c ${bamTag})
  sf=\$(echo "10000000 \$nbreads" | awk '{printf "%.2f", \$1/\$2}')
  
  samtools sort -o ${prefix}_genome1_sorted.bam ${bam1}
  samtools index ${prefix}_genome1_sorted.bam
  bamCoverage -b ${prefix}_genome1_sorted.bam \\
              -o ${prefix}_genome1_rpkm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              --scaleFactor \$sf

  samtools sort -o ${prefix}_genome2_sorted.bam ${bam2}
  samtools index ${prefix}_genome2_sorted.bam
  bamCoverage -b ${prefix}_genome2_sorted.bam \\
              -o ${prefix}_genome2_rpkm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              --scaleFactor \$sf
  """
}


/*****************************************************
 * MultiQC
 */

/*
process get_software_versions {
  output:
  file 'software_versions_mqc.yaml' into software_versions_yaml

  script:
  """
  echo $workflow.manifest.version &> v_main.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  multiqc --version &> v_multiqc.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflow_summary_mqc {
  when:
  !params.skip_multiqc

  output:
  file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/asmapping'
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
    !params.skip_multiqc

    input:
    file splan from ch_splan.collect()
    file metadata from ch_metadata.ifEmpty([])
    file multiqc_config from ch_multiqc_config    
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    file splan
    file "*report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName + "_rnaseq_report" : "--filename rnaseq_report"
    metadata_opts = params.metadata ? "--metadata ${metadata}" : ""
    splan_opts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
    isPE = params.singleEnd ? 0 : 1
    
    modules_list = "-m custom_content -m preseq -m rseqc -m bowtie2 -m hisat2 -m star -m tophat -m cutadapt -m fastqc"
    modules_list = params.counts == 'featureCounts' ? "${modules_list} -m featureCounts" : "${modules_list}"  
    modules_list = params.counts == 'HTseqCounts' ? "${modules_list} -m htseq" : "${modules_list}"  
 
    """
    stats2multiqc.sh ${splan} ${params.aligner} ${isPE}
    ##max_read_nb="\$(awk -F, 'BEGIN{a=0}(\$1>a){a=\$3}END{print a}' mq.stats)"
    median_read_nb="\$(sort -t, -k3,3n mq.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) print(a[x-1]+a[x])/2; else print a[x-1];}')"
    mqc_header.py --name "RNA-seq" --version ${workflow.manifest.version} ${metadata_opts} ${splan_opts} --nbreads \${median_read_nb} > multiqc-config-header.yaml
    multiqc . -f $rtitle $rfilename -c $multiqc_config -c multiqc-config-header.yaml $modules_list
    """    
}
*/


/*
 * Sub-routine
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
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
        log.info "[rnaseq] Pipeline Complete"
    }else{
        log.info "[rnaseq] FAILED: $workflow.runName"
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
