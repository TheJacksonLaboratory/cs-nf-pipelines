process SAMTOOLS_VIEW {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools_view' }", pattern:"*.bam", mode:'copy'

  input:
      tuple val(sampleID), file(sam)
      val(view_string)
      val(filename)

  output:
      tuple val(sampleID), file("*.bam"), emit: bam

  script:
      log.info "----- Samtools View Running on: ${sampleID} -----"

    """
    samtools view ${view_string} ${sam} > ${sampleID}_${filename}.bam 
    """
}


