process CONPAIR_PILEUP {
  tag "$sampleName"

  cpus 1
  memory 4.GB
  time '10:00:00'
  container 'quay.io/jaxcompsci/conpair:v0.2'

  input:
  tuple val(sampleID), val(sampleName), file(bam), file(bai)
  val(type)

  output:
  tuple val(sampleID), val(sampleName), file("*pileup.txt"), emit: pileup

  script:
  """
  python2 /Conpair-0.2/scripts/run_gatk_pileup_for_sample.py -B ${bam} -O ${sampleName}_${type}_pileup.txt --reference ${params.ref_fa} --markers /Conpair-0.2/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed
  """
}

// Marker file inputs: `--markers /Conpair-0.2/data/markers/...` are located in the container. 