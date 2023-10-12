process EMASE_GET_COMMON_ALIGNMENT {
    tag "$sampleID"

    cpus 1
    memory 90.GB
    time 2.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

    input:
    tuple val(sampleID), path(emase_files)

    output:
    tuple val(sampleID), file("*.merged.emase.h5"), emit: common_emase_h5

    script:
    
    emase_list = emase_files.collect { "$it" }.join(' -i ')

    """
    emase get-common-alignments  -i ${emase_list} -o ${sampleID}.merged.emase.h5
    """

    stub:
    """
    touch ${sampleID}.merged.emase.h5
    """
}

/*
 Usage: emase get-common-alignments [OPTIONS]

 get the common alignments

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --emase-file  -i      FILE     EMASE file to compress, can seperate files by "," or have multiple -i [default: None] [required]                                                                                                                                                    │
│    --output      -o      FILE     EMASE file with unique reads [default: None]                                                                                                                                                                                                        │
│    --comp-lib    -c      TEXT     compression library to use [default: zlib]                                                                                                                                                                                                          │
│    --verbose     -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                                         │
│    --help                         Show this message and exit.                                                                                                                                                                                                                         │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/
