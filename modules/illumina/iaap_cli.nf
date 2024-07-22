process IAAP_CLI {
  tag "$idat_folder"
  cpus = 4
  memory { idat_folder.size() < 60.GB ? 8.GB : 24.GB }
  time { idat_folder.size() < 60.GB ? '03:00:00' : '12:00:00' }
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
  
  container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'

  publishDir "${params.pubdir}/${params.organize_by=='idat' ? "$idat_folder" + '/results' : 'iaap_cli'}", pattern:"*.log", mode:'copy'

  input:
  path idat_folder from params.idat_folder
  path output_dir from params.output_dir
  path bpm_file from params.bpm_file
  path egt_file from params.egt_file

  output:
  path "$output_dir/iaap_cli.log", emit: iaap_cli_log

  script:
  """
  mkdir -p $output_dir
  chmod a+w $output_dir

  echo "Running IAAP_CLI with BPM file: $bpm_file and EGT file: $egt_file" > $output_dir/iaap_cli.log

  /usr/local/bin/iaap-cli/iaap-cli gencall \
      $bpm_file \
      $egt_file \
      $output_dir \
      --idat-folder $idat_folder \
      --output-gtc >> $output_dir/iaap_cli.log 2>&1
  """
}

workflow {
  IAAP_CLI()
}