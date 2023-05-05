process GBRS_RECONSTRUCT  {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/mikewlloyd/gbrs_test:latest'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.npz", mode: 'copy', enabled: params.keep_intermediate
    
    input:
    tuple val(sampleID), path(tpm)

    output:
    tuple val(sampleID), file("*genoprobs.npz"), emit: genoprobs_npz
    tuple val(sampleID), file("*genotypes.npz"), emit: genotypes_npz
    tuple val(sampleID), file("*genotypes.tsv"), emit: genotypes_tsv

    script:

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
    touch ${sampleID}.genoprobs.npz
    touch ${sampleID}.genotypes.npz
    touch ${sampleID}.genotypes.tsv
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
