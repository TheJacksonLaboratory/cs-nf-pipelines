process G2GTOOLS_EXTRACT {
    tag "$strain"

    cpus 1
    memory 6.GB
    time '02:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/g2gtools:74926ad'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.fa', mode:'copy'

    input:
    tuple val(strain), path(final_fasta), path(db)
    val(extract_type)

    output:
    tuple val(strain), path("*.fa"), emit: extracted_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    g2gtools extract ${debug_run} -i ${final_fasta} -db ${db} --${extract_type} 2> ${strain}.${params.genome_version}.${extract_type}.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.${extract_type}.fa
    """

}
/* 
NOTE: The above script is hard-coded for extraction of regions from database files. Additional options are available. 
      At this time, only db to fasta is used. If additional functionality is required from this module, additional coding, 
      and variables will be required. See help text below for these options / input formats. Additional conditional logic 
      or sub-modules of the extract function may also be required to avoid parameter conflicts. 
*/ 

/*
 Usage: g2gtools extract [OPTIONS]

 Extract subsequence from a fasta file given a region

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --fasta               -i       FILE     Fasta file to extract from [default: None] [required]                                                                                                                                                                                      │
│    --bed                 -b       FILE     BED file name [default: None]                                                                                                                                                                                                              │
│    --database            -db      FILE     Database file name, use with --genes, --transcripts, --exons [default: None]                                                                                                                                                               │
│    --genes                                 Extract genes from --database                                                                                                                                                                                                              │
│    --transcripts                           Extract transcripts from --database                                                                                                                                                                                                        │
│    --exons                                 Extract exons from --database                                                                                                                                                                                                              │
│    --identifier          -id      TEXT     Fasta identifier [default: None]                                                                                                                                                                                                           │
│    --region              -r       TEXT     Region to extract in chromosome:start-end format [default: None]                                                                                                                                                                           │
│    --vci                 -c       FILE     VCI File to use [default: None]                                                                                                                                                                                                            │
│    --vci-reverse         -R                Reverse the direction of the VCI file                                                                                                                                                                                                      │
│    --complement                            Complement the extracted sequence                                                                                                                                                                                                          │
│    --reverse                               Reverse the extracted sequence                                                                                                                                                                                                             │
│    --reverse-complement                    Reverse-complement the extracted sequence                                                                                                                                                                                                  │
│    --raw                                   Shows just the extracted sequences                                                                                                                                                                                                         │
│    --verbose             -v       INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                                │
│    --help                                  Show this message and exit.                                                                                                                                                                                                                │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
