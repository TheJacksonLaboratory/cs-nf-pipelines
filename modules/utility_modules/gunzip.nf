process GUNZIP {

  tag "$sampleID"

  cpus 1  
  memory { 5.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container "quay.io/jaxcompsci/py3_perl_pylibs:v2"

  input:
  tuple val(sampleID), path(reads)

  output:
  tuple val(sampleID), path("*.{fastq,fq}"), emit: gunzip_fastq
  shell:

  '''
  if [[ !{reads[0]} =~ ".gz" ]];
  then
    gunzip -c !{reads[0]} > !{reads[0].baseName}
  else
    mv !{reads[0]} input_!{reads[0]}
  fi
  if [[ !{reads[1]} =~ ".gz" ]];
  then
    gunzip -c !{reads[1]} > !{reads[1].baseName}
  else
   mv !{reads[1]} input_!{reads[1]}
  fi
  '''
}