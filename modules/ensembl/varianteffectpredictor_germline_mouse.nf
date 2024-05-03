process VEP_GERMLINE {
  tag "$sampleID"

  cpus = 4
  memory = 15.GB
  time = '10:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'ensemblorg/ensembl-vep:release_110.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'vep' }", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(vcf), file(idx)

  output:
  tuple val(sampleID), file("*_vep_annotated.vcf.gz"), emit: vcf_gz
  tuple val(sampleID), file("*_vep_annotated.vcf.gz.tbi"), emit: tbi
  tuple val(sampleID), file("*_vep_annotated.vcf"), emit: vcf

  script:

  """
  vep \
  --input_file ${vcf} \
  --output_file ${sampleID}_germline_vep_annotated.vcf \
  --fork ${task.cpus} \
  --buffer_size 50000 \
  --format vcf \
  --no_stats \
  --no_escape \
  --offline \
  --assembly GRCm39 \
  --species mus_musculus \
  --cache \
  --dir_cache ${params.vep_cache_directory} \
  --refseq \
  --exclude_predicted \
  --fasta ${params.vep_fasta} \
  --symbol \
  --hgvs \
  --check_existing \
  --vcf \
  --pick_allele_gene

  bgzip -c ${sampleID}_germline_vep_annotated.vcf > ${sampleID}_germline_vep_annotated.vcf.gz
  tabix ${sampleID}_germline_vep_annotated.vcf.gz

  """
}

// VEP Cache setup: 

// singularity pull --name vep.sif docker://ensemblorg/ensembl-vep:release_110.1
// singularity exec vep.sif INSTALL.pl -c /PATH_TO_VEP/vep -a cfp -s mus_musculus_refseq -y GRCm39 
