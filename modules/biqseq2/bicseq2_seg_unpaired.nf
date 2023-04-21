process BICSEQ2_SEG_UNPAIRED {
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/bicseq2:v3'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/callers' : 'biqseq2' }", pattern:"{*.txt,*.png}", mode:'copy'

  input:
  tuple val(sampleID), file(individual_tumor_norm_bin_files), val(meta), val(tumor_name)

  output:
  tuple val(sampleID), file("*.bicseq2.png"), val('no_idx'), val(meta), val(params.na12878_sampleName), val(tumor_name), val('bicseq2'), emit: bicseq2_png
  tuple val(sampleID), file("*.bicseq2.txt"), val('no_idx'), val(meta), val(params.na12878_sampleName), val(tumor_name), val('bicseq2'), emit: bicseq2_sv_calls

  script:

  tumor_norm_list = individual_tumor_norm_bin_files.collect { "$it" }.join(' ')

  """

  python3 \
  ${projectDir}/bin/pta/bicseq2_seg_config_writer_unpaired.py \
  --tumor-norms ${tumor_norm_list} \
  --seg-bicseq2-config ${params.bicseq2_chromList} \
  --out-file configuration_file.txt \
  --pair-id ${sampleID}

  perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
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





