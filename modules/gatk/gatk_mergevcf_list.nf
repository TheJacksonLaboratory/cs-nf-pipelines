process GATK_MERGEVCF_LIST {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(list)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    // memory needs to be set explicitly
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" MergeVcfs \
    -R ${params.ref_fa} \
    -I ${list} \
    -O ${sampleID}_GATKcombined_raw.vcf
    """
}
