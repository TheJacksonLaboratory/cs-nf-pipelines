process SNORLAX {
    tag "used sleep. It was super effective."

    cpus = 1
    time = '00:00:30'
    memory = 1.GB
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    // errorStrategy {(task.attempt < 5) ? { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }.call() : 'finish'}
    // errorStrategy 'retry'
    // maxRetries 6

    script:
    """
        echo "START"
        sleep 160
        echo "DONE"
    """
}

// This module is for testing module / nextflow failures only. It is not intended for production level pipeline use. 
