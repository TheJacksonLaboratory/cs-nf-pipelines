process GATK_GENOTYPEGVCF {
    tag "$sampleID"
    
    cpus 4
    array 23
    memory 30.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    // publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy'

    input:
    tuple val(sampleID), path(gvcf), path(idx)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    // memory needs to be set explicitly
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp -XX:ParallelGCThreads=${task.cpus}" GenotypeGVCFs \
    -R ${params.ref_fa} \
    -V ${gvcf} \
    -O ${gvcf.baseName}_genotyped.vcf
    """
}
