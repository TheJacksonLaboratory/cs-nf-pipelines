process DO_TRANSITION_PROBABILITIES {

    cpus 2
    memory 80.GB
    time 2.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v4'

    publishDir "${params.pubdir}", pattern: '*.h5', mode:'copy'
    publishDir "${params.pubdir}", pattern: '*.tsv', mode:'copy'
    publishDir "${params.pubdir}/transprob_matrix_plots", pattern: '*.pdf', mode:'copy'

    output:
    path('*.F.h5'), emit: female_h5_file
    path('*.M.h5'), emit: male_h5_file
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
