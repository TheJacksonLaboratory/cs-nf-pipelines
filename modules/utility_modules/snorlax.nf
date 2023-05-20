process SNORLAX {
  tag "used sleep. It was super effective."

  cpus = 1
  time = '00:00:30'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
  
  script:
  """
    echo "START"
    sleep 160
    echo "DONE"
  """
}
