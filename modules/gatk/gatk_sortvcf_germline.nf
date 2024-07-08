process GATK_SORTVCF_GERMLINE {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '05:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(list)
    val(gvcf)

    output:
    tuple val(sampleID), file("*.vcf"), file("*.idx"), emit: vcf_idx, optional: true

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = list.collect { "-I $it" }.join(' ')

    if (gvcf=='gvcf'){
        output_suffix='g.vcf'
    }
    else{
        output_suffix='vcf'
    }

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" SortVcf  \
    -SD ${params.ref_fa_dict} \
    ${inputs} \
    -O ${sampleID}_merged.${output_suffix}
    """
}
