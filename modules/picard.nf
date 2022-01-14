// part A
process PICARD_ALN_METRICS_A {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.outdir}/picard", pattern: "*.bam", mode: 'copy'
  publishDir "${params.outdir}/picard", pattern: "*.bai", mode: 'copy'

  input:
  tuple val(sampleID), file(read_groups)
  tuple val(sampleID), file(genome_sorted_bam)
  file(picard_dict)
 
  output:
  tuple val(sampleID), file("*group_reorder.bam"), emit: reordered_sorted_bam
  tuple val(sampleID), file("*group_reorder.bai") 
  
  script:
  log.info "----- Picard Alignment Metrics Running on: ${sampleID} -----"

  """
  picard AddOrReplaceReadGroups \
  INPUT=${genome_sorted_bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_groups) \
  CREATE_INDEX=true

  echo "picard 4a.1"

  picard ReorderSam \
  INPUT=${sampleID}_genome_bam_with_read_groups.bam \
  OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
  SEQUENCE_DICTIONARY=${params.picard_dict} \
  CREATE_INDEX=true

  echo "picard 4a.2"

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

  publishDir "${params.outdir}/picard", pattern: "*.txt", mode: 'copy'

  input:
  tuple sampleID, file(reordered_sorted_bam)

  output:
  file "*.*"
  tuple sampleID, file("*metrics.txt"), emit: picard_metrics

  script:
  log.info "----- Alignment Metrics Running on: ${sampleID} -----"

  if (params.read_prep == "stranded" && params.gen_org == "human")

    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /picard.jar SortSam \
    SO=coordinate \
    INPUT=${reordered_sorted_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
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
    INPUT=${reordered_sorted_bam} \
    OUTPUT=${sampleID}_reorder_sort.bam \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

    java -Djava.io.tmpdir=$TMPDIR -Xmx4g -jar /picard.jar CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=NONE \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  else if (params.gen_org == "mouse" && params.reads == "PE")

    """
    bamtools stats -insert -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
    """

  else if (params.gen_org == "mouse" && params.reads != "PE")

    """
    bamtools stats -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
    """

  }
