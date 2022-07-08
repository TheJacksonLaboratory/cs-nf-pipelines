process BAM_COVERAGE_BIGWIG {
  tag "$sampleID"

  cpus 8
  memory 10.GB
  time '04:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'deeptools' }", pattern: "*.bigwig", mode: 'copy'
  container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

  input:
  tuple val(sampleID), file(processed_bams)

  output:
  tuple val(sampleID), file("*.bigwig")

  script:
  log.info "----- Running deeptools bamCoverage bigwig on ${sampleID} -----"
  """
  bamCoverage \
  --numberOfProcessors $task.cpus \
  --binSize 10 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize ${params.effective_genome_size} \
  --bam ${processed_bams[0]} \
  --outFileFormat bigwig \
  --outFileName ${sampleID}.bigwig
  """
}
