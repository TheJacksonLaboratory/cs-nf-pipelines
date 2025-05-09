process GATK_HAPLOTYPECALLER {
    tag "$sampleID"

    cpus = 1
    memory = 30.GB
    time = '15:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.*vcf", mode:'copy'

    input:
    tuple val(sampleID), file(bam), file(bai)
    val(gvcf)

    output:
    tuple val(sampleID), file("*.*vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    if (gvcf=='gvcf'){
    delta="-ERC GVCF"
    output_suffix='gvcf'
    }
    else{
    delta="--dbsnp ${params.dbSNP} -stand-call-conf ${params.call_val}"
    output_suffix='vcf'
    }


    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" HaplotypeCaller  \
    -R ${params.ref_fa} \
    -I ${bam} \
    -O ${sampleID}_haplotypecaller_raw.${output_suffix} \
    -L ${params.target_gatk} \
    -ploidy ${params.ploidy_val} \
    ${delta} \
    --max-reads-per-alignment-start 0 \
    --recover-all-dangling-branches
    """
}
