process GATK_GETPILEUPSUMMARIES {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(bam), path(bai), val(sample_name)

    output:
    tuple val(sampleID), val(meta), path("*pileups.table"), emit: pileup_summary

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" GetPileupSummaries \
    -I ${bam} \
    -V ${params.contam_ref} \
    -L ${params.contam_ref} \
    -O ${sampleID}.pileups.table
    """
}
