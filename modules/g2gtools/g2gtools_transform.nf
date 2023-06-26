process G2GTOOLS_TRANSFORM {
    tag "$strain"

    cpus 8
    memory 25.GB
    time '02:30:00'

    container 'quay.io/jaxcompsci/g2gtools:74926ad'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.fa', mode:'copy'

    input:
    tuple val(strain), path(patched_fasta), path(vci), path(tbi)

    output:
    tuple val(strain), path("*.fa"), emit: final_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    g2gtools transform -p ${task.cpus} ${params.region} ${params.bed} ${debug_run} -i ${patched_fasta} -c ${vci} -o ${strain}.${params.genome_version}.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.fa
    """

}

/*
 Usage: g2gtools transform [OPTIONS]

 Incorporate indels onto the input sequence

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
