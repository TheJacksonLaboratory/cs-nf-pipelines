process ANNOTATE_BICSEQ2_CNV {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.bed", mode: 'copy'

    input:
        //BICSEQ2_SEG.out.bicseq2_sv_calls
        tuple val(sampleID), file(bicseq2_calls), val(no_idx), val(meta), val(normal_name), val(tumor_name), val(bicseq2)
        val(chrom_list)

    output:
        tuple val(sampleID), file("${sampleID}_cnv_annotated_final.bed"), val(normal_name), val(tumor_name), emit: bicseq_annot
        tuple val(sampleID), file("${sampleID}_cnv_annotated_supplemental.bed"), val(normal_name), val(tumor_name), emit: bicseq_annot_suppl

    script:
        listOfChroms = chrom_list.collect { "$it" }.join(',')

        """
        Rscript ${projectDir}/bin/pta/annotate-cnv.r \
            --cnv=${bicseq2_calls} \
            --caller="bicseq2" \
            --tumor=${tumor_name} \
            --normal=${normal_name} \
            --cytoband=${params.cytoband} \
            --db_names="DGV,1000G,COSMIC" \
            --db_files=${params.dgv},${params.thousandG},${params.cosmicUniqueBed} \
            --cancer_census=${params.cancerCensusBed} \
            --ensembl=${params.ensemblUniqueBed} \
            --allowed_chr=${listOfChroms} \
            --overlap_fraction=0.8 \
            --out_file_main=${sampleID}_cnv_annotated_final.bed \
            --out_file_supplemental=${sampleID}_cnv_annotated_supplemental.bed

        """
}
