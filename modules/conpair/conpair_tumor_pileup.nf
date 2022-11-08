process CONPAIR_TUMOR_PILEUP {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '10:00:00'
  container 'quay.io/jaxcompsci/conpair:v0.2'

  input:
  tuple val(sampleID), file(tumor_bam), file(normal_bai), val(meta)

  output:
  tuple val(sampleID), file("*_tumor_pileup.txt"), val(meta), emit: tumor_pileup

  script:
  """
  python2 /Conpair-0.2/scripts/run_gatk_pileup_for_sample.py -B ${tumor_bam} -O ${sampleID}_tumor_pileup.txt --reference ${params.ref_fa} --markers /Conpair-0.2/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed
  """
}
