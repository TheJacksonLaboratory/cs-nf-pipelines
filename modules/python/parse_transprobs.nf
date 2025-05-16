process PARSE_TRANSITION_PROBABILITIES {
    tag "${generation} - ${sex}"
    
    cpus 2
    memory 2.GB
    time 1.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'hdfgroup/h5py:2.7.0'

    publishDir "${params.pubdir}", pattern: '*.npz', mode:'copy'

    input:
    tuple path(h5), val(generation)
    val(sex)

    output:
    path('*.npz'), emit: npz_files

    script:

    """
    python ${projectDir}/bin/gbrs/parse_h5_transprob_to_npz.py -t ${h5} -s ${sex} -g ${generation} -l "${params.haplotype_list}"
    """

    stub:
    """
    touch "tranprob.DO.G1.F.npz"
    """
}
