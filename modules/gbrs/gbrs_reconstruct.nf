process GBRS_RECONSTRUCT  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*", mode: 'copy', enabled: "${ suffix == 'merged' || params.read_type == 'SE' ? true : false }"

    input:
    tuple val(sampleID), path(tpm)

    output:
    tuple val(sampleID), file("*"), emit: compressed_emase_h5

    script:

    // output_name = suffix == 'merged' ? "${sampleID}.merged.compressed.emase.h5" : "${bam[0].baseName}.compressed.emase.h5"

    """

    cp ${params.base_ref_index_fai} ref.fa.fai

    gbrs reconstruct \
        -e ${tpm} \
        -t ${params.trans_prob_dir}/tranprob.DO.G${params.sample_generation}.${params.sample_sex}.npz \
        -x ${params.emission_prob_avecs} \
        -g ${params.gene_position_file} \
        -o ${sampleID} \
        -c ${params.gbrs_expression_threshold} \
        -s ${params.gbrs_sigma} 
    """

    stub:
    """
    touch ${sampleID}
    """
}

/*
usage: gbrs reconstruct [-h] -e EXPRFILE -t TPROBFILE [-x AVECFILE]
                        [-g GPOSFILE] [-c EXPR_THRESHOLD] [-s SIGMA]
                        [-o OUTBASE]

optional arguments:
  -h, --help         show this help message and exit
  -e EXPRFILE
  -t TPROBFILE
  -x AVECFILE
  -g GPOSFILE
  -c EXPR_THRESHOLD
  -s SIGMA
  -o OUTBASE

 sigma defaults to 0.12.
 expr_threshold defaults to 1.5. 

*/
