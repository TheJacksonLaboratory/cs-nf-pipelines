process MAKE_CUSTOM_TRANSCRIPTOME {
    tag "$custom_gene_fasta"

    cpus 1
    memory 15.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    publishDir "${params.pubdir}/custom_transcriptome", mode:'copy'

    input:
    tuple path(ref_fa), path(ref_gtf), path(custom_gene_fasta)

    output:
    path("*.fasta"), emit: concat_fasta
    path("${custom_gene_fasta.baseName}.gtf"), emit: custom_gtf
    path("${ref_fa.baseName}_${custom_gene_fasta.baseName}.gtf"), emit: concat_gtf

    script:
        """
        cat ${ref_fa} ${custom_gene_fasta} > ${ref_fa.baseName}_${custom_gene_fasta.baseName}.fasta
        python ${projectDir}/bin/rnaseq/fasta_to_gtf.py -i ${custom_gene_fasta} -o ${custom_gene_fasta.baseName}.gtf
        cat ${ref_gtf} ${custom_gene_fasta.baseName}.gtf > ${ref_fa.baseName}_${custom_gene_fasta.baseName}.gtf
        """
}
