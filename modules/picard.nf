//Picard Tools
// Will need to depricate pamA and pamB for new approach

process PICARD_SORTSAM{
  tag "sampleID"

  // Slurm Options
  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  // Container
  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // Publish Directory
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*_aln.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(bwa_mem)

  output:
  tuple val(sampleID), file(*_aln.ba*), emit: picard_sortsam

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"

  """
  picard SortSam \
  SO=coordinate \
  INPUT=${sam} \
  OUTPUT=${sampleID}_aln_sortsam.bam  \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true
  """
}

process PICARD_MARKDUPLICATES{
  tag "sampleID"

  // Slurm Options
  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  // Container
  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // Publish Directory
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*_aln.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(sorted_sam)

  output:
  tuple val(sampleID), file("*_dedup.ba*"), emit: bam_dedup
  tuple val(sampleID), file("*metrics.dat"), emit: picard_metrics

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"

  """
  picard MarkDuplicates \
  I=${sorted_sam ? *.bam} \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.dat \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT
  """
}

// part A
process PICARD_ALN_METRICS_A {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(read_groups)
  tuple val(sampleID), file(genome_sorted_bam)


  output:
  tuple val(sampleID), file("*group_reorder.bam"), emit: reordered_sorted_bam
  tuple val(sampleID), file("*group_reorder.bai"), emit: reordered_sorted_bai

  script:
  log.info "----- Picard Alignment Metrics Running on: ${sampleID} -----"

  """
  picard AddOrReplaceReadGroups \
  INPUT=${genome_sorted_bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_groups) \
  CREATE_INDEX=true

  picard ReorderSam \
  INPUT=${sampleID}_genome_bam_with_read_groups.bam \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  SEQUENCE_DICTIONARY=${params.picard_dict} \
  CREATE_INDEX=true
  """

  }

// part B
process PICARD_ALN_METRICS_B {

  // human only for mouse see bamtools
  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(reordered_sorted_bam)

  output:
  tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics

  script:
  log.info "----- Alignment Metrics B Human Running on: ${sampleID} -----"

  if (params.read_prep == "stranded")

    """
    picard SortSam \
    SO=coordinate \
    INPUT=${reordered_sorted_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    picard CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  else if (params.read_prep != "stranded")

    """
    picard SortSam \
    SO=coordinate \
    INPUT=${reordered_sorted_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    picard CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=NONE \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  }
