process CALC_MTDNA_FILTER_CHRM {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*_mtDNA_Content.txt", mode: 'copy'
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(rmdup_bam_file)
  tuple val(sampleID), file(rmdup_bai_file)

  output:
  tuple val(sampleID), file("*.sorted.rmDup.rmChrM.bam"), emit: rmChrM_bam
  tuple val(sampleID), file("*.sorted.rmDup.rmChrM.bam.bai"), emit: rmChrM_bai
  tuple val(sampleID), file("*_mtDNA_Content.txt"), emit: mtdna_log

  shell:
  log.info "----- Calculate %mtDNA and Filter Mitochondrial Reads on ${sampleID} -----"
  '''
  mtReads=$(samtools idxstats !{rmdup_bam_file} | grep 'MT' | cut -f 3)
  totalReads=$(samtools idxstats !{rmdup_bam_file} | awk '{SUM += $3} END {print SUM}')

  if [ $mtReads >0 ]
  then
    mtReads=$(echo $mtReads)
  else
    mtReads=$(echo 0)
  fi

  echo 'mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' >> !{sampleID}_mtDNA_Content.txt

  samtools view -@ !{params.threads} -h !{rmdup_bam_file} \
  | grep -v MT \
  | samtools sort -@ !{params.threads} -O bam \
  -o !{sampleID}.sorted.rmDup.rmChrM.bam \
  && samtools index !{sampleID}.sorted.rmDup.rmChrM.bam
  '''

}
