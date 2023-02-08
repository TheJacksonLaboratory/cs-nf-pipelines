process STRELKA2 {
  tag "$meta.patient"

  cpus = 4
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/strelka2:v2.9.3'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'strelka' }", pattern:"*.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), file(candidateSmallIndels), file(candidateSmallIndels_tbi)

  output:
  tuple val(sampleID), file("*indels.vcf.gz"), file("*indels.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), val('strelka'), emit: strelka_indel_vcf_tbi
  tuple val(sampleID), file("*snvs.vcf.gz"), file("*snvs.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), val('strelka'), emit: strelka_snv_vcf_tbi

  script:

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  # configure strelka
  python /usr/local/bin/configureStrelkaSomaticWorkflow.py \
  --normalBam ${normal_bam} \
  --tumorBam ${tumor_bam} \
  --callRegions ${params.callRegions} \
  --referenceFasta ${params.ref_fa} \
  --indelCandidates ${candidateSmallIndels} \
  --config ${params.strelka_config} \
  --runDir ${sampleID}

  # execute strelka
  python ${sampleID}/runWorkflow.py \
  --mode local \
  --job ${task.cpus} \
  --memGb ${my_mem}

  mv ${sampleID}/results/variants/somatic.snvs.vcf.gz ${sampleID}_strelka_somatic.snvs.vcf.gz
  mv ${sampleID}/results/variants/somatic.snvs.vcf.gz.tbi ${sampleID}_strelka_somatic.snvs.vcf.gz.tbi
  mv ${sampleID}/results/variants/somatic.indels.vcf.gz ${sampleID}_strelka_somatic.indels.vcf.gz
  mv ${sampleID}/results/variants/somatic.indels.vcf.gz.tbi ${sampleID}_strelka_somatic.indels.vcf.gz.tbi

  """
}