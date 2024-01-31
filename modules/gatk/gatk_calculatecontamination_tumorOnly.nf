process GATK_CALCULATECONTAMINATION {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(pileup_table)

    output:
    tuple val(sampleID), path("*contamination.table"), path("*segments.txt"), emit: contam_segments

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus}" CalculateContamination \
    -I ${pileup_table} \
    --tumor-segmentation ${sampleID}.segments.txt \
    -O ${sampleID}.contamination.table
    """
}
