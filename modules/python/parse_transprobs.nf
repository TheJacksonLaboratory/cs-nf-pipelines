process PARSE_TRANSITION_PROBABILITIES {

    cpus 2
    memory 5.GB
    time 3.hour

    container 'hdfgroup/h5py:2.7.0'

    publishDir "${params.pubdir}", pattern: '*.npz', mode:'copy'

    input:
    tuple path(female_h5), path(male_h5)

    output:
    path('*.npz'), emit: npz_files

    script:
    """
    python ${projectDir}/bin/gbrs/parse_h5_transprob_to_npz.py -m ${male_h5} -f ${female_h5} -g ${params.num_generations}
    """

    stub:
    """
    touch "tranprob.DO.G1.F.npz"
    touch "tranprob.DO.G1.M.npz"
    """

}
