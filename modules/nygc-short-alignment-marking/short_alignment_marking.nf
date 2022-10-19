process SHORT_ALIGNMENT_MARKING {
  tag "$sampleID"

  cpus 2
  memory 4.GB
  time '10:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'short_alignment_marking' }", pattern:"*.marked.bam", mode:'copy'
  
  input:
  tuple val(sampleID), file(aligned_bam)

  output:
  tuple val(sampleID), file("*.marked.bam"), emit: marked_bam

  script:
  log.info "----- Short Alignment Marking Running on: ${sampleID} -----"
  // parses the bam file and marks as unmapped a read with alignment length below a user-defined threshold. Reads are not filtered from the bam file but kept as unmapped.
  """
  /projects/compsci/jgeorge/harshpreet_CS/nextflow/pipeline/ngs-ops/nygc_pipeline/tools/nygc-short-alignment-marking/filter_bam -I ${alignmed_bam} -A1 30 -A2 30 -o ${sampleID}.marked.bam

  """

}

