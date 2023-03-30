process GATK_MUTECT2 {
  tag "$sampleID"

  cpus = 4
  memory = 15.GB
  time {15.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'broadinstitute/gatk:4.4.0.0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), file("*_somatic.vcf.gz"), file("*_somatic.vcf.gz.tbi"), file("*.stats"), emit: vcf_tbi_stats

  script:
  //Estimate somatic variants using Mutect2
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus}" Mutect2 \
    -R ${params.ref_fa} \
    -I ${tumor_bam} \
    --germline-resource ${params.gnomad_ref} \
    --panel-of-normals ${params.pon_ref} \
    --genotype-germline-sites true \
    --genotype-pon-sites true \
    --pileup-detection \
    --dont-use-soft-clipped-bases false \
    -L ${params.target_gatk} \
    --native-pair-hmm-threads 4 \
    --annotation QualByDepth \
    --annotation RMSMappingQuality \
    --annotation FisherStrand \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest \
    --min-base-quality-score 20 \
    -O ${sampleID}_mutect2_somatic.vcf.gz
  """
}

/*
As of v4.1, there is no longer a need to specify the tumor sample name with -tumor. You need only specify the normal sample name with -normal, if you include a normal.

Starting with v4.0.4.0, GATK recommends the default setting of --af-of-alleles-not-in-resource, which the tool dynamically adjusts for different modes. 
tumor-only calling sets the default to 5e-8, tumor-normal calling sets it to 1e-6 and mitochondrial mode sets it to 4e-3. 
For previous versions, the default was 0.001, the average heterozygosity of humans. 
For other organisms, change --af-of-alleles-not-in-resource to 1/(ploidy*samples in resource).
*/