process GATK_MUTECT2 {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), path(interval), val(interval_index)

    output:
    tuple val(sampleID), path("*_somatic.vcf.gz"), val(meta), val(normal_name), val(tumor_name), val('mutect2'), emit: vcf
    tuple val(sampleID), path("*_somatic.vcf.gz.tbi"), emit: tbi
    tuple val(sampleID), path("*.stats"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" Mutect2 \
    -R ${params.ref_fa} \
    -I ${tumor_bam} \
    -tumor ${tumor_name} \
    -I ${normal_bam} \
    -normal ${normal_name} \
    -L ${interval} \
    --native-pair-hmm-threads 4 \
    -O ${sampleID}_${interval_index}_somatic.vcf.gz
    """
}
