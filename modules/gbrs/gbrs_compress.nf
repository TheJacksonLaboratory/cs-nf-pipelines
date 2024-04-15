process GBRS_COMPRESS {
    tag "$sampleID"

    cpus 1
    memory { params.read_type == 'SE' ? 12.GB : 250.GB }
    time 5.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:v1.0.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/emase' : 'emase' }", pattern: "*.h5", mode: 'copy'

    input:
    tuple val(sampleID), path(h5)
    val(suffix)

    output:
    tuple val(sampleID), file("*.compressed.emase.h5"), emit: compressed_emase_h5

    script:
    output_name = suffix == 'merged' ? "${sampleID}.merged.compressed.emase.h5" : "${sampleID}.compressed.emase.h5"

    """
    gbrs compress -i ${h5} -o ${output_name}
    """

    stub:
    """
    touch ${output_name}
    """
}

/*
 Usage: gbrs compress [OPTIONS]

 compress EMASE format alignment incidence matrix

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --emase-file  -i      FILE     EMASE file to compress, can seperate files by "," or have multiple -i [default: None] [required]                                                                                                                                 │
│ *  --output      -o      FILE     name of the compressed EMASE file [default: None] [required]                                                                                                                                                                     │
│    --comp-lib    -c      TEXT     compression library to use [default: zlib]                                                                                                                                                                                       │
│    --verbose     -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                      │
│    --help                         Show this message and exit.                                                                                                                                                                                                      │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/
