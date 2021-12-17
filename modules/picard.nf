// part A
process PICARD_ALN_METRICS_A {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'java_samtools_python_R_picard_bamtools.sif'

  input:
  tuple sampleID, file(gen_sort_bam), file(read_grp)

  output:
  tuple sampleID, file("*group_reorder.bam")

  script:
  log.info "----- Picard Alignment Metrics Running on: ${sampleID} -----"

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar AddOrReplaceReadGroups \
  INPUT=${gen_sort_bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_group.bam \
  SORT_ORDER=coordinate \
  $(cat $read_grp) \
  CREATE_INDEX=true

  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar ReorderSam \
  INPUT=${sampleID}_genome_bam_with_read_group.bam \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  REFERENCE=${params.ref_fa} \
  CREATE_INDEX=true
  """
  }

// part B
process PICARD_ALN_METRICS_B {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'java_samtools_python_R_picard_bamtools.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(reord_sort_bam)

  output:
  file "*.*"
  tuple sampleID, file("*metrics.txt")

  script:
  log.info "----- Alignment Metrics Running on: ${sampleID} -----"

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
