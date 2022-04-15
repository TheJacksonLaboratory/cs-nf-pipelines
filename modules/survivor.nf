process SURVIVOR{
  tag "$sampleID"

  cpus = 1
  time = '2:00:00'
  memory = '100.GB'
  clusterOptions = '-q batch'
  maxRetries = 1
  errorStrategy = 'retry'

  container ''
  input:
  val(vcf_paths)
		// uses following params
    // params.names
		// params.surv_dist
		// params.surv_supp
		// params.surv_type
		// params.surv_strand
		// params.surv_min
	output:
	tuple val(sampleID), file("*.vcf")

  script:
  if('illumina'){.BDLM}
  if('pacbio'){.PS}
	"""
	SURVIVOR merge ${vcf_paths} ${params.surv_dist} ${surv_supp} ${surv_type} ${surv_strand} 0 ${surv_min} ${sample_name}_mergedCall_$process.vcf
	"""
}
