process G2GTOOLS_GTF2DB {
    tag "$strain"

    cpus 1
    memory 1.GB
    time '02:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/g2gtools:74926ad'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.gtf.db', mode:'copy'

    input:
    tuple val(strain), path(gtf)

    output:
    tuple val(strain), path("*.gtf.db"), emit: db

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    g2gtools gtf2db ${debug_run} -i ${gtf} -o ${strain}.${params.genome_version}.gtf.db
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.gtf.db
    """

}

/*
 Usage: g2gtools fasta-format [OPTIONS]

 Convert a GTF file to a G2G DB file

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input    -i      FILE     GTF file [default: None] [required]                                                                                                                                                                                                                    │
│    --output   -o      FILE     Name of output file [default: None]                                                                                                                                                                                                                    │
│    --verbose  -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                                            │
│    --help                      Show this message and exit.                                                                                                                                                                                                                            │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
