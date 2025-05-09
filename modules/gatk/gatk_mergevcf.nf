process GATK_MERGEVCF {
    tag "$sampleID"

    cpus 1
    memory 48.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.workflow == 'amplicon' || params.workflow == 'amplicon_generic' ? params.keep_intermediate : true

    input:
    tuple val(sampleID), file(snp_vcf), file(indel_vcf)
    val(suffix)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" MergeVcfs \
    -R ${params.ref_fa} \
    -I ${snp_vcf} \
    -I ${indel_vcf} \
    -O ${sampleID}_${suffix}.vcf
    """
}
