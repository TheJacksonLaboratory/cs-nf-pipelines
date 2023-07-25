process GBRS_INTERPOLATE  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-547132f'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.npz", mode: 'copy'

    input:
    tuple val(sampleID), path(genoprobs_npz)

    output:
    tuple val(sampleID), file("*gbrs.interpolated.genoprobs.npz"), emit: interpolated_genoprobs


    script:

    """
    cp ${params.base_ref_index_fai} ref.fa.fai

    gbrs interpolate \
        -i ${genoprobs_npz} \
        -g ${params.genotype_grid} \
        -p ${params.gene_position_file} \
        -o ${sampleID}.gbrs.interpolated.genoprobs.npz
    """

    stub:
    """
    touch ${sampleID}.gbrs.interpolated.genoprobs.npz
    """
}

/*
 Usage: gbrs interpolate [OPTIONS]

 interpolate probability on a decently-spaced grid

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --genoprob-file  -i      FILE     genotype probability file [default: None] [required]                                                                                                                                                                          │
│    --grid-file      -g      FILE     grid file (i.e, ref.genome_grid.64k.txt) [default: None]                                                                                                                                                                      │
│    --gpos-file      -p      FILE     meta information for genes (chrom, id, location) [default: None]                                                                                                                                                              │
│    --output         -o      FILE     output file in GBRS quant format [default: None]                                                                                                                                                                              │
│    --verbose        -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                   │
│    --help                            Show this message and exit.                                                                                                                                                                                                   │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
