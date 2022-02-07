//Picard Tools

process PICARD_SORTSAM {
  tag "sampleID"

  // Slurm Options
  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  // Container
  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // Publish Directory
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*_sortsam.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(sam)

  output:
  tuple val(sampleID), file("*_sortsam.bam"), emit: bam
  tuple val(sampleID), file("*_sortsam.bai"), emit: bai

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"

  """
  picard SortSam \
  SO=coordinate \
  INPUT=${sam} \
  OUTPUT=${sampleID}_sortsam.bam  \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true
  """
}

process PICARD_MARKDUPLICATES {
  tag "sampleID"

  // Slurm Options
  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  // Container
  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // Publish Directory
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.dat", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_dedup.bam"), emit: dedup_bam
  tuple val(sampleID), file("*_dedup.bai"), emit: dedup_bai
  tuple val(sampleID), file("*.dat"), emit: dedup_metrics

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"

  """
  picard MarkDuplicates \
  I=${bam} \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.dat \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT
  """
}


process PICARD_COLLECTHSMETRICS {
  tag "sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.*", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*Metrics.txt"), emit: hsmetrics

  script:
  log.info "----- Picard CollectHsMetrics Running on: ${sampleID} -----"

  """
  picard CollectHsMetrics \
  INPUT=${bam} \
  OUTPUT=${sampleID}_CoverageMetrics.txt \
  BAIT_INTERVALS=${params.bait_picard} \
  TARGET_INTERVALS=${params.target_picard} \
  REFERENCE_SEQUENCE=${params.ref_fa} \
  VALIDATION_STRINGENCY=SILENT
  """
}

process PICARD_ADDORREPLACEREADGROUPS {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(read_groups)
  tuple val(sampleID), file(bam)


  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- Picard Add or Replace Read Groups Running on: ${sampleID} -----"

  """
  picard AddOrReplaceReadGroups \
  INPUT=${bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_groups) \
  CREATE_INDEX=true
  """

  }
  process PICARD_REORDERSAM {

    tag "sampleID"

    cpus 1
    memory 8.GB
    time '12:00:00'
    clusterOptions '-q batch'

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'

    input:
    tuple val(sampleID), file(bam)


    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.bai"), emit: bai

    script:
    log.info "----- Picard Alignment Metrics Running on: ${sampleID} -----"

    """
    picard ReorderSam \
    INPUT=${bam} \
    OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
    SEQUENCE_DICTIONARY=${params.picard_dict} \
    CREATE_INDEX=true
    """

    }

process PICARD_COLLECTRNASEQMETRICS {

  // human only for mouse see bamtools
  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.pdf", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics
  tuple val(sampleID), file("*.pdf")

  script:
  log.info "----- Alignment Metrics B Human Running on: ${sampleID} -----"

  if (params.read_prep == "stranded")

    """
    picard CollectRnaSeqMetrics \
    I=${bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  else if (params.read_prep != "stranded")

    """
    picard CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=NONE \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  }
