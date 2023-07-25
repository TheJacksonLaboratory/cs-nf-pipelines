process G2GTOOLS_PATCH {
    tag "$strain"

    cpus 8
    memory 8.GB
    time '02:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/g2gtools:74926ad'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.patched.fa', mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(strain), path(vci), path(tbi)

    output:
    tuple val(strain), path("*.patched.fa"), emit: patched_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    g2gtools patch -p ${task.cpus} ${params.region} ${params.bed} ${debug_run} -i ${params.primary_reference_fasta} -c ${vci} -o ${strain}.${params.genome_version}.patched.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.patched.fa
    """

}

/*
 Usage: g2gtools patch [OPTIONS]

 Patch SNPs onto the reference sequence

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input    -i      FILE     Fasta file to extract from [default: None] [required]                                                                                                                                                                                                  │
│ *  --vci      -c      FILE     VCI File to use [default: None] [required]                                                                                                                                                                                                             │
│    --output   -o      FILE     Name of output file [default: None]                                                                                                                                                                                                                    │
│    --bed      -b      FILE     BED file name [default: None]                                                                                                                                                                                                                          │
│    --region   -r      TEXT     Region to extract in chromosome:start-end format [default: None]                                                                                                                                                                                       │
│    --reverse                   Reverse the direction of VCI file                                                                                                                                                                                                                      │
│    --bgzip                     compress and index output                                                                                                                                                                                                                              │
│    --verbose  -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                                            │
│    --help                      Show this message and exit.                                                                                                                                                                                                                            │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
