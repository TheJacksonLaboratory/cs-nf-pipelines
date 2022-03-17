process LUMPY_MAP {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*_lumpy_aligned_bwa.bam"), emit: aligned
  tuple val(sampleID), file("*_lumpy_discordant.bam"), emit: discordant

  script:
  log.info "----- Lumpy Mapping Running on: ${sampleID} -----"

  """
  # Clipped_rc reads mapping to Genome
  samtools sort -n ${bam} -o ${sampleID}_lumpy.bam -@ ${task.cpus}

  # manual read group info
  samtools view -h ${sampleID}_lumpy.bam \
  | samblaster --acceptDupMarks --excludeDups --addMateTags \
  --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
  | samtools view -@ ${task.cpus} -S -b - > ${sampleID}_lumpy_aligned_bwa.bam

  # Extract the discordant pairedend alignments
  samtools view -@ ${task.cpus} -b -F 1294 ${sampleID}_aligned_lumpy.bam > ${sampleID}_lumpy_discordant.bam
  """
}
process LUMPY_EXTRACT_SPLITS {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- Lumpy Extract Running on: ${sampleID} -----"

  """
  samtools view -h ${bam} \
  | extractSplitReads_BwaMem -i stdin \
  | samtools view -Sb - > ${sampleID}_splitters_unsorted.bam
  """
}
process LUMPY_CALL_SV {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

  input:
  tuple val(sampleID), file(aligned_bam)
  tuple val(sampleID), path(discordant_bam)
  tuple val(sampleID), path(split_bam)


  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  // move to config
  params.exclude_regions = "/ref_data/mm10.gaps.centro_telo.scafold.exclude.bed"

  """
  RG_ID=$(samtools view -H ${aligned_bam} | grep '^@RG' | sed "s/.*ID:\\([^\\t]*\\).*/\\1/g")
  samtools view -r "${RG_ID}" ${aligned_bam} | tail -n+100000 > pre_metrics 2>/dev/null
  metrics=$(cat pre_metrics | ${params.pairend_distro} -r 150 -X 4 -N 10000 -o ${sampleID}.histo) 2>/dev/null \
    && [ $? = 141 ] && echo 'metrics to pairend_distro had exitcode: '$?;
  mean=$(echo "${metrics}" | cut -d " " -f 1)
  mean=$(echo "${mean}"    | cut -d ":" -f 2)
  std_dev=$(echo "${metrics}" | cut -d " " -f 2)
  std_dev=$(echo "${std_dev}" | cut -d ":" -f 2)

  lumpy \
  -mw 4 \
  -x ${params.exclude_regions} \
  -pe id:"${RG_ID}",bam_file:${discordant_bam},histo_file:${sampleID}.histo,mean:"${mean}",stdev:"${std_dev}",read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
  -sr id:"${RG_ID}",bam_file:${split_bam},back_distance:10,weight:1,min_mapping_threshold:20 \
  > ${lumpy_vcf}

  vcfSort.sh ${lumpy_vcf} ${lumpy_sort_vcf}
  """
}
