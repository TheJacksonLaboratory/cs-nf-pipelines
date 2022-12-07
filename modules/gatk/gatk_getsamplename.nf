process GATK_GETSAMPLENAME {
  tag "$meta.sample"

  cpus = 1
  memory = 1.GB
  time = '00:05:00'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), val(meta), file(bam), file(bai)
  
  output:
  tuple val(sampleID), stdout, emit: sample_name
 
  script:
  """
  gatk GetSampleName \
    -I ${bam} \
    -O sample_name.txt

  cat sample_name.txt
  """
}