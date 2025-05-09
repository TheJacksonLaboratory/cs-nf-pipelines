process GOOGLE_DEEPVARIANT {

  cpus 2
  memory 50.GB
  time {10.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  container 'docker://google/deepvariant:latest'

  input:
  tuple val(sampleID), file(bam), file(bai), val(chrom), val(ind), val(sex)

  output:
  tuple val(sampleID), file("*.vcf.gz"), file("*.gvcf.gz"), file("*.vcf.gz.tbi"), file("*.gvcf.gz.tbi"), emit: vcf_channel

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
