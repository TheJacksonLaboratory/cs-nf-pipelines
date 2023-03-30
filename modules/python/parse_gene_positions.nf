process PARSE_GENE_POSITONS {

    cpus 2
    memory 5.GB
    time 10.minutes 

    container 'hdfgroup/h5py:2.7.0'

    publishDir "${params.pubdir}", pattern: '*.npz', mode:'copy'

    input:
    path(gene_pos_tsv)

    output:
    path('*.npz'), emit: npz_file

    script:
    """
    python ${projectDir}/bin/gbrs/parse_ref_gene_pos_file.py -i ${gene_pos_tsv} -t ${params.emase_gene2transcript} -n "ref.gene_pos.ordered_ensBuild_${params.ensembl_build}"
    """

    stub:
    """
    touch "ref.gene_pos.ordered.npz"
    """

}


