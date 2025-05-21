process GATK_GENOTYPE_GVCF {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.*vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(vcf_index), path(interval), val(index)

    output:
    tuple val(sampleID), file("*.*vcf"), file("*.idx"), path(interval), val(index), emit: vcf_idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" GenotypeGVCFs  \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -O ${sampleID}_${index}_genotypedGVCFs.vcf \
    -L ${interval}
    """
}
