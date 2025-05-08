process BICSEQ2_SEG_UNPAIRED {
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/bicseq2:v3'
  publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern:"{*.txt,*.png}", mode:'copy'

  input:
  tuple val(sampleID), file(individual_tumor_norm_bin_files), val(meta), val(tumor_name)

  output:
  tuple val(sampleID), file("*.bicseq2.png"), val('no_idx'), val(meta), val(params.na12878_sampleName), val(tumor_name), val('bicseq2'), emit: bicseq2_png
  tuple val(sampleID), file("*.bicseq2.txt"), val('no_idx'), val(meta), val(params.na12878_sampleName), val(tumor_name), val('bicseq2'), emit: bicseq2_sv_calls

  script:

  tumor_norm_list = individual_tumor_norm_bin_files.collect { "$it" }.join(' ')

  scale = params.bicseq2_no_scaling ? "--noscale" : ""

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
  ${scale} \
  configuration_file.txt \
  ${sampleID}.bicseq2.txt
  
  """

  stub:
  """
  wget -O ${sampleID}.bicseq2.png https://raw.githubusercontent.com/TheJacksonLaboratory/cs-nf-test/main/pta/human/fizzbang--t_bang--n_fizz.bicseq2.png
  wget -O ${sampleID}.bicseq2.txt https://raw.githubusercontent.com/TheJacksonLaboratory/cs-nf-test/main/pta/human/fizzbang--t_bang--n_fizz.bicseq2.txt
  """
}
// Stub run is to test lower coverage sample data.
