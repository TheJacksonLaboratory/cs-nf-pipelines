process GBRS_INTERPOLATE  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/mikewlloyd/gbrs_test:latest'

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
usage: gbrs interpolate [-h] -i PROBFILE [-g GRIDFILE] [-p GPOSFILE]
                        [-o OUTFILE]

optional arguments:
  -h, --help   show this help message and exit
  -i PROBFILE
  -g GRIDFILE
  -p GPOSFILE
  -o OUTFILE
*/
