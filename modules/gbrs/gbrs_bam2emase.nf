process GBRS_BAM2EMASE {
    tag "$sampleID"

    cpus 1
    memory { params.read_type == 'SE' ? 12.GB : 45.GB }
    time 6.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gbrs' }", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    gbrs bam2emase -i ${bam} \
                -m ${params.transcripts_info} \
                -h ${params.gbrs_strain_list} \
                -o ${bam.baseName}.emase.h5
    """

    stub:
    """
    touch ${bam.baseName}.emase.h5
    """
}


/*
 Usage: gbrs bam2emase [OPTIONS]

 convert a BAM file to EMASE format

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --alignment-file  -i      FILE     bam file to convert [default: None] [required]                                                                                                                                                                               │
│ *  --haplotype-char  -h      TEXT     haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C [default: None] [required]                                                                                                                  │
│ *  --locus-ids       -m      FILE     filename for the locus (usually transcripts) info [default: None] [required]                                                                                                                                                 │
│    --output          -o      FILE     EMASE file (hdf5 format) [default: None]                                                                                                                                                                                     │
│    --delim           -d      TEXT     delimiter string between locus and haplotype in BAM file [default: _]                                                                                                                                                        │
│    --verbose         -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                  │
│    --help                             Show this message and exit.                                                                                                                                                                                                  │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/