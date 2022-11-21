process GATKv3_5_HAPLOTYPECALLER {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'broadinstitute/gatk3:3.5-0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'gatk' }", pattern: "*.gvcf", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)

  output:
  tuple val(meta), file("*.gvcf"), emit: normal_germline_vcf
  tuple val(meta), file("*.gvcf.idx"), emit: normal_germline_index

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
  -T HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${normal_bam} \
  -o ${sampleID}_variants_raw.gvcf \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  -ERC GVCF
  """
}