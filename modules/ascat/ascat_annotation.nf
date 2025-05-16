process ASCAT_ANNOTATION {
    tag "$sampleID"

    cpus = 1
    memory = 24.GB
    time = '01:30:00'
    errorStrategy = { (task.exitStatus == 140) ? { log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish' }.call() : 'finish' }

    container 'quay.io/jaxcompsci/ascat:v3.1.3'

    publishDir "${params.pubdir}/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(segments_raw), path(ploidy)

    output:
    tuple val(sampleID), val(meta), path("*.segments_raw.extend.txt"), emit: seg_extended
    tuple val(sampleID), val(meta), path("*.ensgene_cnvbreak.txt"), emit: ensembl_annot
    tuple val(sampleID), val(meta), path("*.png"), emit: png

    script:
    gender = meta.gender == 'XX' ? 'female' : 'male'
    """
    perl ${projectDir}/bin/cnv_array/segment_raw_extend.pl ${segments_raw} ${ploidy} ${params.chrArm} ${gender}
    perl ${projectDir}/bin/cnv_array/annotate_ensembl_genes.pl ${sampleID}.segments_raw.extend.txt ${params.cnvGeneFile}
    R CMD BATCH --slave "--args ${sampleID}.segments_raw.extend.txt ${sampleID} ./ " ${projectDir}/bin/cnv_array/seg_plot.R
    """
}
