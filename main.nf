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
      --reads 'READS'               Path to input data (must be surrounded with quotes)
      --samplePlan 'SAMPLEPLAN'     Path to sample plan input file (cannot be used with --reads)
      --genome 'BUILD'              Name of genome reference
      -profile PROFILE              Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

    Sequencing:
      --singleEnd                   Specifies that the input is single end reads

    Genotype:
      --maternal
      --paternal

    Mapping:
      --nmask
      --aligner 'MAPPER'            Tool for read alignments ['star', 'hisat2', 'tophat2']. Default: 'star'

    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index 'PATH'           Path to STAR index
      --hisat2_index 'PATH'         Path to HiSAT2 index
      --tophat2_index 'PATH'        Path to TopHat2 index
      --saveAlignedIntermediates    Save the intermediate files from the Aligment step  - not done by default

    Other options:
      --metadata 'FILE'             Add metadata file for multiQC report
      --outdir 'PATH'               The output directory where the results will be saved
      -w/--work-dir 'PATH'          The temporary directory where intermediate data will be saved
      --email 'MAIL'                Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name 'NAME'                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

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

Channel.fromPath("${params.genome}")
       .ifEmpty { exit 1, "Reference Genome not found: ${params.genome}" }
       .set { fastaGenome }

if ( params.metadata ){
   Channel
       .fromPath( params.metadata )
       .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
       .set { ch_metadata }
}

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; raw_reads_strandness; save_strandness}
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; raw_reads_strandness; save_strandness}
   }
   params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; raw_reads_strandness; save_strandness}
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; raw_reads_strandness; save_strandness}
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; raw_reads_strandness; save_strandness}
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
  if(params.star_index) summary['STAR Index'] = params.star_index
} else if(params.aligner == 'tophat2') {
  summary['Aligner'] = "Tophat2"
  if(params.bowtie2_index) summary['Tophat2 Index'] = params.bowtie2_index
} else if(params.aligner == 'hisat2') {
  summary['Aligner'] = "HISAT2"
  if(params.hisat2_index) summary['HISAT2 Index'] = params.hisat2_index
}
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
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
 * Prepare Genome
 */


process prepareReferenceGenome {

  input:
  file reference from fastaGenome.collect()

  output:
  file "*.fa" into referenceGenome

  if (params.nmask){
  """
  SNPsplit_genome_preparation --strain ${params.paternal} --strain2 ${params.maternal} --reference_genome ${REF_DIR} --vcf_file ${VCF} --nmasking
  """
  }else{
  script:
  """
  SNPsplit_genome_preparation --strain ${params.paternal} --strain2 ${params.maternal} --reference_genome ${REF_DIR} --vcf_file ${VCF} --no_nmasking
  """
  }
}

/*
 * MultiQC
 */

process get_software_versions {
  output:
  file 'software_versions_mqc.yaml' into software_versions_yaml

  script:
  """
  echo $workflow.manifest.version &> v_main.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  fastqc --version &> v_fastqc.txt
  STAR --version &> v_star.txt
  tophat2 --version &> v_tophat2.txt
  hisat2 --version &> v_hisat2.txt
  preseq &> v_preseq.txt
  infer_experiment.py --version &> v_rseqc.txt
  read_duplication.py --version &> v_read_duplication.txt
  featureCounts -v &> v_featurecounts.txt
  htseq-count -h | grep version  &> v_htseq.txt
  picard MarkDuplicates --version &> v_markduplicates.txt || true
  samtools --version &> v_samtools.txt
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
  section_href: 'https://gitlab.curie.fr/rnaseq'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
  """.stripIndent()
}

process multiqc {
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

    report_fields['skipped_poor_alignment'] = skipped_poor_alignment

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

    if(skipped_poor_alignment.size() > 0){
        log.info "[rnaseq] WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor alignment scores!"
    }

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
