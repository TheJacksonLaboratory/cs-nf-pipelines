process SEQUENZA_ANNOTATE {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '00:45:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/perl-vcftools-vcf:0.1.16--pl5321hdfd78af_4'

    publishDir "${params.pubdir}/${sampleID + '/sequenza_cnv'}", pattern:"*.txt", mode:'copy'

    input:
    tuple val(sampleID), val(meta), path(segment_file)

    output:
    tuple val(sampleID), path("*.txt"), emit: annotated_segment

    script:
    """
    perl ${projectDir}/bin/wes/ensembl_annotation.pl ${segment_file} ${params.ensembl_database}
    """
}
