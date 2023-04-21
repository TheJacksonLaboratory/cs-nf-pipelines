process BICSEQ2_SEG {
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/bicseq2:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'biqseq2' }", pattern:"{*.txt,*.png}", mode:'copy'

  input:
  tuple val(sampleID), file(individual_normal_norm_bin_files), file(individual_tumor_norm_bin_files), val(meta), val(normal_name), val(tumor_name)

  output:
  tuple val(sampleID), file("*.bicseq2.png"), val('no_idx'), val(meta), val(normal_name), val(tumor_name), val('bicseq2'), emit: bicseq2_png
  tuple val(sampleID), file("*.bicseq2.txt"), val('no_idx'), val(meta), val(normal_name), val(tumor_name), val('bicseq2'), emit: bicseq2_sv_calls

  script:

  normal_norm_list = individual_normal_norm_bin_files.collect { "$it" }.join(' ')
  tumor_norm_list = individual_tumor_norm_bin_files.collect { "$it" }.join(' ')


  """

  python3 \
  ${projectDir}/bin/pta/bicseq2_seg_config_writer.py \
  --normal-norms ${normal_norm_list} \
  --tumor-norms ${tumor_norm_list} \
  --seg-bicseq2-config ${params.bicseq2_chromList} \
  --out-file configuration_file.txt \
  --pair-id ${sampleID}

  perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
  --control \
  --tmp ${sampleID} \
  --fig ${sampleID}.bicseq2.png \
  --title ${sampleID} \
  --lambda 4 \
  configuration_file.txt \
  ${sampleID}.bicseq2.txt
  
  """

  stub:
  """
  touch ${sampleID}.bicseq2.png
  touch ${sampleID}.bicseq2.txt
  """
}





