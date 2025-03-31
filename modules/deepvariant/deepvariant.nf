process GOOGLE_DEEPVARIANT {

  cpus 1
  memory 50.GB
  time '18:00:00'

  container 'docker://google/deepvariant:latest'

  input:
  tuple val(sampleID), file(bam), file(bai), val(chrom), val(ind), val(sex)

  output:
  tuple val(sampleID), val(ind), val(sex), file("*_dv.vcf"), emit: vcf
  tuple val(sampleID), val(ind), val(sex), file("*_dv.g.vcf"), emit: gvcf

  script:
  
  """
  if [ ${sex} == "M" ]; then
    run_deepvariant \
    --model_type WGS \
    --ref ${params.ref_fa} \
    --haploid_contigs X,Y \
    --reads ${bam} \
    --output_vcf ${sampleID}_${chrom}.vcf \
    --sample_name ${sampleID} \
    --num_shards 4 \
    --regions ${chrom} \
    --intermediate_results_dir ${params.pubdir}/${sampleID}/${chrom} \
    --output_gvcf ${sampleID}_${chrom}.g.vcf \
    --verbosity 1
  fi

  if [[ ${sex} == "F" || ${sex} == "U" ]]; then 
    run_deepvariant \
    --model_type WGS \
    --ref ${params.ref_fa} \
    --reads ${bam} \
    --output_vcf ${sampleID}_${chrom}_dv.vcf \
    --sample_name ${sampleID} \
    --num_shards 4 \
    --regions ${chrom} \
    --intermediate_results_dir ${params.pubdir}/${sampleID}/${chrom} \
    --output_gvcf ${sampleID}_${chrom}_dv.g.vcf \
    --verbosity 1
  fi
  
  """
}