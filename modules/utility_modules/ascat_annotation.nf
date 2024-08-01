process ASCAT_ANNOTATION {

    tag "$sampleID"

    cpus = 1
    memory = 24.GB
    time = '01:30:00'
    errorStrategy = { (task.exitStatus == 140) ? { log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish' }.call() : 'finish' }

    container 'quay.io/biocontainers/ascat:3.1.1--r43hdfd78af_1'
    publishDir "${params.pubdir}/${params.organize_by == 'sample' ? sampleID : 'ascat_annotation'}", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(segments_raw), path(ploidy)

    output:
    tuple val(sampleID), val(meta), path("${sampleID}.segments_raw.extend.txt"), path("${sampleID}.*"), emit: ascat_annotated

    script:
    """
    perl \${projectDir}/bin/cnv_array/${sampleID}.segment_raw_extend.pl ${segments_raw} ${ploidy} ${params.chrArm} ${meta}
    perl \${projectDir}/bin/cnv_array/annotate_ensembl_genes.pl ${sampleID}.segments_raw.extend.txt ${params.cnvGeneFile}
    R CMD BATCH --slave "--args ${sampleID}.segments_raw.extend.txt ${sampleID} ./ " \${projectDir}/seg_plot.R
    """
}
