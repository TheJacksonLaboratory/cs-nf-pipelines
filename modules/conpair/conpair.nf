process CONPAIR {
  tag "$meta.patient"

  cpus 1
  memory 4.GB
  time '10:00:00'
  container 'quay.io/jaxcompsci/conpair:v0.2'
  errorStrategy 'ignore'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'conpair' }", pattern:"*.txt", mode:'copy'
  
  input:
  tuple val(meta), file(tumor_pileup), file(normal_pileup)

  output:
  tuple val(meta), file("*_concordance.txt"), emit: concordance
  tuple val(meta), file("*_contamination.txt"), emit: contamination

  script:

  //Estimate concordance and contamination estimator for tumor–normal pairs
  //Verifying concordance between two samples (tumor and normal)
  //Estimating contamination level in both the tumor and the normal:
  """
  python2 /Conpair-0.2/scripts/verify_concordance.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${meta.patient}_concordance.txt

  python2 /Conpair-0.2/scripts/estimate_tumor_normal_contamination.py -T ${tumor_pileup} -N ${normal_pileup} --outfile ${meta.patient}_contamination.txt
  """

  stub:
  """
  touch ${meta.patient}_concordance.txt
  touch ${meta.patient}_contamination.txt
  """
}
