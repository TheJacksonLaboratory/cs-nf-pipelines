process MANTA {
  tag "$meta.patient"

  cpus = 4
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/manta:v1.5.0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'manta' }", pattern:"*.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), file("*candidateSmallIndels.vcf.gz"), file("*candidateSmallIndels.vcf.gz.tbi"), emit: manta_smallindel_vcf_tbi
  tuple val(sampleID), file("*diploidSV.vcf.gz"), file("*diploidSV.vcf.gz.tbi"), emit: manta_diploidsv_tbi
  tuple val(sampleID), file("*somaticSV.vcf.gz"), file("*somaticSV.vcf.gz.tbi"), emit: manta_somaticsv_tbi
  tuple val(sampleID), file("*candidateSV.vcf.gz"), file("*candidateSV.vcf.gz.tbi"), emit: manta_candidatesv_tbi


  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  # configure manta
  python configManta.py \
  --normalBam ${normal_bam} \
  --tumorBam ${tumor_bam} \
  --referenceFasta ${params.ref_fa} \
  --callRegions ${params.callRegions} \
  --runDir ${sampleID}

  # execute manta
  python ${sampleID}/runWorkflow.py -j ${task.cpus} \
  --mode local \
  --memGb ${my_mem}
  
  mv ${sampleID}/results/variants/candidateSmallIndels.vcf.gz ${sampleID}_manta_candidateSmallIndels.vcf.gz
  mv ${sampleID}/results/variants/candidateSmallIndels.vcf.gz.tbi ${sampleID}_manta_candidateSmallIndels.vcf.gz.tbi
  mv ${sampleID}/results/variants/diploidSV.vcf.gz ${sampleID}_manta_diploidSV.vcf.gz
  mv ${sampleID}/results/variants/diploidSV.vcf.gz.tbi ${sampleID}_manta_diploidSV.vcf.gz.tbi
  mv ${sampleID}/results/variants/somaticSV.vcf.gz ${sampleID}_manta_somaticSV.vcf.gz
  mv ${sampleID}/results/variants/somaticSV.vcf.gz.tbi ${sampleID}_manta_somaticSV.vcf.gz.tbi
  mv ${sampleID}/results/variants/candidateSV.vcf.gz ${sampleID}_manta_candidateSV.vcf.gz
  mv ${sampleID}/results/variants/candidateSV.vcf.gz.tbi ${sampleID}_manta_candidateSV.vcf.gz.tbi
  """
}