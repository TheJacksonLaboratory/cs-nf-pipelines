process GBRS_PLOT  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-16c7011'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.pdf", mode: 'copy'

    input:
    tuple val(sampleID), path(interpolated_genoprobs)

    output:
    tuple val(sampleID), file("*gbrs.plotted.genome.pdf"), emit: genotype_pdf


    script:

    """
    cp ${params.base_ref_index_fai} ref.fa.fai
    cp ${params.founder_hex_colors} founder.hexcolor.info

    gbrs plot \
        -i ${interpolated_genoprobs} \
        -o ${sampleID}.gbrs.plotted.genome.pdf \
        -n ${sampleID} 
    """

    stub:
    """
    touch ${sampleID}.gbrs.plotted.genome.pdf
    """
}

/*
 Usage: gbrs plot [OPTIONS]

 plot a reconstructed genome

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --genoprob-file  -i      FILE     EMASE genoprobs file [default: None] [required]                                                                                                                                                                               │
│    --output         -o      FILE     name of output file [default: None]                                                                                                                                                                                           │
│    --format         -f      TEXT     output file format [default: pdf]                                                                                                                                                                                             │
│    --sample_name    -n      TEXT     name of the sample [default: None]                                                                                                                                                                                            │
│    --verbose        -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                   │
│    --help                            Show this message and exit.                                                                                                                                                                                                   │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

*/
