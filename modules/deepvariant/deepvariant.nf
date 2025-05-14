process GOOGLE_DEEPVARIANT {

  cpus 2
  memory 50.GB
  time {10.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  container 'docker://google/deepvariant:1.9.0'

  input:
  tuple val(sampleID), file(bam), file(bai), val(chrom), val(sex)

  output:
  tuple val(sampleID), path("*.vcf.gz"), path("*.gvcf.gz"), path("*.vcf.gz.tbi"), path("*.gvcf.gz.tbi"), emit: vcf_channel

  script:
  
  """
  if [ ${sex} == "M" ]; then
    run_deepvariant \
    --model_type WGS \
    --ref ${params.ref_fa} \
    --haploid_contigs X,Y \
    --reads ${bam} \
    --output_vcf ${sampleID}_${chrom}.vcf.gz \
    --sample_name ${sampleID} \
    --num_shards 4 \
    --regions ${chrom} \
    --output_gvcf ${sampleID}_${chrom}.gvcf.gz \
    --verbosity 1 
  fi

  if [[ ${sex} == "F" || ${sex} == "U" ]]; then 
    run_deepvariant \
    --model_type WGS \
    --ref ${params.ref_fa} \
    --reads ${bam} \
    --output_vcf ${sampleID}_${chrom}.vcf.gz \
    --sample_name ${sampleID} \
    --num_shards 4 \
    --regions ${chrom} \
    --output_gvcf ${sampleID}_${chrom}.gvcf.gz \
    --verbosity 1
  fi
  
  """
}
