process DO_TRANSITION_PROBABILITIES {

    cpus 2
    memory 40.GB
    time 48.hour

    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}", pattern: '*.h5', mode:'copy'
    publishDir "${params.pubdir}", pattern: '*.tsv', mode:'copy'
    publishDir "${params.pubdir}/transprob_matrix_plots", pattern: '*.pdf', mode:'copy'

    output:
    tuple path('*.F.h5'), path('*.M.h5'), emit: h5_files
    path('*.tsv'), emit: gene_list_tsv
    path('*.pdf'), emit: pdf_files

    script:
    """
    Rscript ${projectDir}/bin/gbrs/gene_bp_to_cM_to_transprob.R --ensembl_build ${params.ensembl_build} --num_generation ${params.num_generations} --output_prefix "tranprob.genes.DO."
    """

    stub:
    """
    touch "tranprob.genes.DO.G0-1.F.h5"
    touch "tranprob.genes.DO.G0-1.M.h5"
    touch "tranprob.genes.DO.G0-1.F.genPlots.pdf"
    touch "tranprob.genes.DO.G0-1.M.genPlots.pdf"
    """

}
