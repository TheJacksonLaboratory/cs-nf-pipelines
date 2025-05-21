process ANNOTATE_GENES_SV {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    input:
        tuple val(sampleID), file(annot_sv_bedpe), val(normal_name), val(tumor_name)
        val(suppl_switch)

    output:
        tuple val(sampleID), file("*.manta_gridss_sv_annotated_genes*.bed"), val(normal_name), val(tumor_name), emit: annot_sv_genes_bedpe

    script:

    if (suppl_switch == "main")
        """
        Rscript ${projectDir}/bin/pta/annotate-bedpe-with-genes.r \
            --ensembl=${params.ensemblUniqueBed} \
            --cancer_census=${params.cancerCensusBed} \
            --bedpe=${annot_sv_bedpe} \
            --out_file=${sampleID}.manta_gridss_sv_annotated_genes.bed
        """
    else if (suppl_switch == "supplemental")
        """
        Rscript ${projectDir}/bin/pta/annotate-bedpe-with-genes.r \
            --ensembl=${params.ensemblUniqueBed} \
            --cancer_census=${params.cancerCensusBed} \
            --bedpe=${annot_sv_bedpe} \
            --out_file=${sampleID}.manta_gridss_sv_annotated_genes_supplemental.bed \
            --supplemental
        """
}
