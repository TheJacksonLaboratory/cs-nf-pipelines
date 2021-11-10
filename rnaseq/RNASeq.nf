#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.fastqInputs          = null // absolute path to input fastq(s)
    params.outdir               = null // directory where output will be stored
    params.tmpdir               = null // directory where files will be stored while pipeline is running
    params.min_pct_hq_reads     = null // minimum percentage of high quality reads; currently set as 50%                                                     
    params.reads                = null // type of sequencing reads; options are PE (paired end) or SE (single end)                
    params.read_prep            = null // type of RNA preparation protocol; options stranded or non_stranded                          
    params.seed_length          = null // seed length for RSEM alignment step; currently set as 25
    params.gen_org              = null // species; options are human or mouse                                                           
    params.aligner              = null // aligner used by RSEM; currently set as bowtie2                              
    params.rsem_ref_prefix      = null // absolute path to, and prefix of, reference genome files
    params.ref_fa               = null // absolute path to reference genome .fa file          
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.config run path/to/pipeline.nf -profile slurm, singularity
            (The params.config file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --fastqInputs             absolute path to input fastq(s).
        --outdir                  directory where output will be stored.
        --tmpdir                  directory where temporary files will be held.
        --min_pct_hq_reads        minimum percentage of high quality reads
        --reads                   type of sequencing reads
        --read_prep               type of RNA preparation protocol
        --seed_length             seed length for RSEM alignment step
        --gen_org                 species; options are human or mouse
        --aligner                 aligner used by RSEM
        --rsem_ref_prefix         absolute path to, and prefix of, reference genome files
        --ref_fa                  absolute path to reference genome .fa file

    Optional:
        -name                   Name for the pipeline run. If not specified Nextflow will
                                automatically generate a random mnemonic.

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.fastqInputs ) {
    exit 1, "Parameter ERROR: fastqInputs ($params.fastqInputs) must be the absolute path to fastq file(s) with _R1 and _R2 fastq filenames."
}
if ( ! params.outdir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}
if ( ! file(params.outdir).exists() ) {
    exit 1, "Parameter ERROR: Missing output directory ($params.outdir) is not found: check if path is correct."
}
if ( ! params.tmpdir ) {
    exit 1, "Parameter ERROR: Temporary directory parameter must be specified."
}
if ( ! file(params.tmpdir).exists() ) {
    exit 1, "Parameter ERROR: Missing temporary directory ($params.tmpdir) is not found: check if path is correct."
}
if ( ! params.ref_fa ) {
    exit 1, "Parameter ERROR: absolute path to reference genome .fa file must be specified."
}
if ( ! params.rsem_ref_prefix ) {
    exit 1, "Parameter ERROR: absolute path to, and prefix of, reference genome files must be specified."
}
if ( ! params.min_pct_hq_reads ) {
    exit 1, "Parameter ERROR: minimum percentage of high quality reads must be specified."
}
if ( ! params.reads ) {
    exit 1, "Parameter ERROR: type of sequencing reads (PE or SE)  must be specified."
}
if ( ! params.read_prep ) {
    exit 1, "Parameter ERROR: type of RNA preparation protocol (stranded or non_stranded) must be specified."
}
if ( ! params.seed_length ) {
    exit 1, "Parameter ERROR: seed length for RSEM alignment step must be specified."
}
if ( ! params.gen_org ) {
    exit 1, "Parameter ERROR: species (human or mouse) must be specified."
}
if ( ! params.aligner ) {
    exit 1, "Parameter ERROR: aligner used by RSEM must be specified."
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
println "Pipeline         = ${workflow.manifest.name}"
println "Description      = ${workflow.manifest.description}"
if(workflow.revision) {
    println "Pipeline Release = ${workflow.revision}"
}
println "Run Name         = ${workflow.runName}"
println "User             = ${workflow.userName}"
println "Config Profile   = ${workflow.profile}"
println "Config Files     = ${workflow.configFiles}"
println "Command Line     = ${workflow.commandLine}"
println "Nextflow Info    = v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
println "Launch dir       = ${workflow.launchDir}"
println "Working dir      = ${workflow.workDir}"
println "Workflow dir     = ${workflow.projectDir}"

// Pipeline Params:
println "Parameters......"
println ".  Output dir                              = ${params.outdir}"
println ".  Temporary dir                           = ${params.tmpdir}"
println ".  FASTQ inputs                            = ${params.fastqInputs}"
println ".  Reference genome                        = ${params.ref_fa}"
println ".  Reference genome files prefix           = ${params.rsem_ref_prefix}"
println ".  Minimum percentage high quality reads   = ${params.min_pct_hq_reads}"
println ".  Type of sequencing reads                = ${params.reads}"
println ".  Type of RNA prep protocol               = ${params.read_prep}"
println ".  Seed length for RSEM alignment step     = ${params.seed_length}"
println ".  Species                                 = ${params.gen_org}"
println ".  Aligner used by RSEM                    = ${params.aligner}"
println "Run Start Time                             = ${workflow.start}"


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def timestamp = new Date().format("yyyyMMdd-HHmmss")
def tmpdir    = params.tmpdir


def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/

fqin = params.fastqInputs.tokenize(",")
fqR1 = file(fqin[0])


// Need to specify absolute path to fastq file(s) because relative paths cause symlink breakage
fqR1p = fqR1.toAbsolutePath().toString()
fqR2p = file(fqin[1]).toAbsolutePath().toString()

def sampleID = ( fqR1.name - fqPairRE )


//~~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1p, fqR2p )
       .toList()
       .set { sample_fastqs_ch }


//~~~~~~~~~~~~~~~~ Primary publishDir == sample_tmpdir ~~~~~
def sample_tmpdir = "${params.tmpdir}/${timestamp}_${sampleID}"


// Step 1: Qual_Stat
process qual_stat {
  tag "sampleID"
  label 'long_mem'
  label 'python_2_7_3'

  publishDir "${sample_tmpdir}_tmp", pattern: "*fastq.gz_stat", mode: 'copy'

  input:
  tuple sampleID, read1, read2 from sample_fastqs_ch

  output:
  tuple sampleID, file("${sampleID}_R{1,2}*filtered_trimmed") into trimmed_fastq
  file "*.fastq.gz_stat"
  tuple sampleID, file("*.fastq.gz_stat") into fq_stats, dummy_fq_stats

  script:
  log.info "-----Qual_Stat running on ${sampleID}-----"

  if (params.reads == "SE"){
    mode_HQ="-S -M"
    inputfq="${read1}"
  }
  if (params.reads == "PE"){
    mode_HQ="-M"
    inputfq="${read1} ${read2}"
  }
  """

  python ${params.filter_trim} $mode_HQ ${params.min_pct_hq_reads}  $inputfq

  """
}


// Step 2: RSEM
process rsem_alignment_exp {
  tag "sampleID"
  label 'long_high_mem'
  label 'rsem'

  publishDir "${sample_tmpdir}_tmp", pattern: "*stats", mode: 'copy'
  publishDir "${sample_tmpdir}_tmp", pattern: "*results*", mode: 'copy'

  input:
  tuple sampleID, file(trimmed) from trimmed_fastq

  output:
  file "*stats"
  file "*results*"
  tuple sampleID, file("*aln.stats") into rsem_stats, dummy_rsem_stats
  tuple sampleID, file("*genome.bam") into genome_sorted_bam
  tuple sampleID, file("*genes.results") into results_genes, dummy_results_genes
  tuple sampleID, file("*isoforms.results") into results_isoforms, dummy_results_isoforms

  script:
  log.info "-----Genome alignment running on ${sampleID}-----"

  if (params.read_prep == "stranded"){
    prob="--forward-prob 0"
  }
  if (params.read_prep == "non_stranded"){
    prob="--forward-prob 0.5"
  }

  if (params.reads == "PE"){
    frag=""
    trimmedfq="--paired-end ${trimmed[0]} ${trimmed[1]}"
  }
  if (params.reads == "SE"){
    frag="--fragment-length-mean 280 --fragment-length-sd 50"
    trimmedfq="${trimmed[0]}"
  }



  """

  rsem-calculate-expression -p 12 --phred33-quals $frag --seed-length ${params.seed_length} $prob --time --output-genome-bam ${params.aligner} $trimmedfq  ${params.rsem_ref_prefix} ${sampleID} 2> ${sampleID}_rsem_aln.stats

  """
  }


Channel.of( sampleID, fqR1p, fqR2p )
    .toList()
    .set { sample_fastqs_ch2 }


// Step 3: Get Read Group Information
process read_group {
  tag "sampleID"
  label 'vshort_mem'
  label 'python_2_7_3'

  publishDir "${sample_tmpdir}_tmp", pattern: "*read_group.txt", mode: 'copy'

  input:
  tuple sampleID, file(read) from sample_fastqs_ch2

  output:
  tuple sampleID, file("*.txt") into read_grp1, read_grp2

  script:
  log.info "-----Read group information determination running on ${sampleID}-----"

  """

  python ${params.read_grp_det} -p -o ${sampleID}_read_group.txt ${read[0]}

  """
  }


// Step 4a: Picard Alignment Metrics, part 1
process picard_aln_metrics_a {
  tag "sampleID"
  label 'med_mem'
  label 'picard_metrics'

  input:
  tuple sampleID, file(gen_sort_bam), file(read_grp) from genome_sorted_bam.join(read_grp1)

  output:
  tuple sampleID, file("*group_reorder.bam") into reordered_sorted_bam1, reordered_sorted_bam2

  script:
  log.info "-----Picard alignment metrics running on ${sampleID}-----"


  """

  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar AddOrReplaceReadGroups \
  INPUT=${gen_sort_bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_group.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_grp) \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar ReorderSam \
  INPUT=${sampleID}_genome_bam_with_read_group.bam \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  REFERENCE=${params.ref_fa} \
  CREATE_INDEX=true

  """
  }

// Step 4b: Picard Alignment Metrics, part 2
process picard_aln_metrics_b {
  tag "sampleID"
  label 'med_mem'
  label 'picard_metrics'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(reord_sort_bam) from reordered_sorted_bam1

  output:
  file "*.*"
  tuple sampleID, file("*metrics.txt") into picard_metrics, dummy_picard_metrics_hsa, dummy_picard_metrics_mmu

  script:
  log.info "-----Alignment metrics running on ${sampleID}-----"

  if (params.read_prep == "stranded" && params.gen_org == "human")

    """

    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar SortSam \
    SO=coordinate \
    INPUT=${reord_sort_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
    I=${reord_sort_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf

    """

  else if (params.read_prep != "stranded" && params.gen_org == "human")

    """

    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar SortSam \
    SO=coordinate \
    INPUT=${reord_sort_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
    I=${reord_sort_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=NONE \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf

    """

  else if (params.gen_org == "mouse" && params.reads == "PE")

    """

    bamtools stats -insert -in ${reord_sort_bam} > ${sampleID}_aln_metrics.txt

    """

  else if (params.gen_org == "mouse" && params.reads != "PE")

    """

    bamtools stats -in ${reord_sort_bam} > ${sampleID}_aln_metrics.txt

    """

  }


// Step 5: Summary Stats (Human samples only)
process summ_stats {
    tag "sampleID"
    label 'xshort_mem'
    label 'perl_R'

    publishDir "${sample_tmpdir}_tmp", pattern: "*stats.txt", mode: 'copy'

    input:
    tuple sampleID, file(mets_stat) from picard_metrics
    tuple sampleID, file(aln_stat) from rsem_stats
    tuple sampleID, file(fq_stat) from fq_stats

    output:
    tuple sampleID, file("*.txt") into summ_stats, dummy_summ_stats

    when:
    params.gen_org == "human"

    script:
    log.info "-----Human Summary Metrics running on ${sampleID}-----"

    if (params.reads == "PE")

      """

      perl ${params.summary_mets_PE} ${fq_stat} ${aln_stat} ${mets_stat} > ${sampleID}_summary_stats.txt

      """

    else if (params.reads == "SE")

      """
      perl ${params.summary_mets_SE} ${fq_stat} ${aln_stat} ${mets_stat}  > ${sampleID}_summary_stats.txt

      """

    }



// Step 6a: GATK Coverage Stats (Human samples only)
process gatk_cov_stats_a {
  tag "sampleID"
  label 'long_mem'
  label 'gatk'


  input:
  tuple sampleID, file(reord_sorted_bam) from reordered_sorted_bam2

  output:
  tuple sampleID, file("*gatk_temp3*") into gatk_stat3
  tuple sampleID, file("*gatk_temp6*") into gatk_stat6

  when:
  params.gen_org == "human"

  script:
  log.info "-----Human GATK coverage stats, part 1 running on ${sampleID}-----"

  """

  samtools index ${reord_sorted_bam}


  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${reord_sorted_bam} \
  -L  ${params.probes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS


  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp4.txt \
  -I ${reord_sorted_bam} \
  -L ${params.ctp_genes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS


  ${params.gatk_form} ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${params.probes}

  ${params.gatk_form} ${sampleID}_gatk_temp4.txt ${sampleID}_gatk_temp5.txt ${sampleID}_gatk_temp6.txt ${params.ctp_genes}

  """
  }



// Step 6b: GATK Coverage Stats (Human samples only)

process gatk_stats_b {
  tag "sampleID"
  label 'long_mem'
  label 'python_2_7_3'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(gatk_temp3), file(gatk_temp6) from gatk_stat3.join(gatk_stat6)


  output:
  file "*exome_interval_avg_median_coverage.bed"
  file "*CCP_interval_avg_median_coverage.bed"
  tuple sampleID, file("*CP_interval_avg_median_coverage.bed") into dummy_avg_med_cov_bed

  when:
  params.gen_org == "human"

  script:
  log.info "-----Human GATK coverage stats, part 2 running on ${sampleID}-----"

  """

  python ${params.cov_calc} ${sampleID}_gatk_temp3.txt ${sampleID}_exome_interval_avg_median_coverage.bed

  python ${params.cov_calc} ${sampleID}_gatk_temp6.txt ${sampleID}_CCP_interval_avg_median_coverage.bed

  """
  }


// Step 7: Move final files to sample directories (human samples)
process transfer_files_hsa {
  tag "sampleID"
  label 'trans_mem'

  input:
  tuple sampleID, file(fastq_stat) from dummy_fq_stats
  tuple sampleID, file(rsem_stat) from dummy_rsem_stats
  tuple sampleID, file(gene_results) from dummy_results_genes
  tuple sampleID, file(isoform_results) from dummy_results_isoforms
  tuple sampleID, file(picard_metrics) from dummy_picard_metrics_hsa
  tuple sampleID, file(sum_stats) from dummy_summ_stats
  tuple sampleID, file(avg_med_coverage) from dummy_avg_med_cov_bed

  when:
  params.gen_org == "human"

  script:
  log.info "-----Moving files to output directory for ${sampleID}-----"

  """

  echo "${sample_tmpdir}" > sampletmpdir.txt

  fullpath=`cat sampletmpdir.txt`


  finfullpath=\$(basename \$fullpath)


  mv "\${fullpath}_tmp" "${params.outdir}/\${finfullpath}/"


  """
  }

// Step 7: Move final files to sample directories (mouse samples)
process transfer_files_mmu {
  tag "sampleID"
  label 'trans_mem'

  input:
  tuple sampleID, file(fastq_stat) from dummy_fq_stats
  tuple sampleID, file(rsem_stat) from dummy_rsem_stats
  tuple sampleID, file(gene_results) from dummy_results_genes
  tuple sampleID, file(isoform_results) from dummy_results_isoforms
  tuple sampleID, file(picard_metrics) from dummy_picard_metrics_mmu

  when:
  params.gen_org == "mouse"

  script:
  log.info "-----Moving files to output directory for ${sampleID}-----"

  """

  echo "${sample_tmpdir}" > sampletmpdir.txt

  fullpath=`cat sampletmpdir.txt`


  finfullpath=\$(basename \$fullpath)


  mv "\${fullpath}_tmp" "${params.outdir}/\${finfullpath}/"


  """
  }

// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Summary ~ ~ ~ ~ ~
workflow.onComplete {
  wfEnd = [:]
  wfEnd['Completed at'] = ${workflow.complete}
  wfEnd['Duration']     = ${workflow.duration}
  wfEnd['Exit status']  = ${workflow.exitStatus}
  wfEnd['Success']      = ${workflow.success}
    if(!workflow.success){
      wfEnd['!!Execution failed'] = ''
      wfEnd['.    Error']   = ${workflow.errorMessage}
      wfEnd['.    Report']  = ${workflow.errorReport}
    }
  Summary.show(wfEnd)
}

