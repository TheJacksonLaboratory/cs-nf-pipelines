process GATKv3_5_HAPLOTYPECALLER {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk3:3.5-0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.gvcf", mode:'copy'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai)

    output:
    tuple val(sampleID), path("*.gvcf"), emit: normal_germline_gvcf
    tuple val(sampleID), path("*.gvcf.idx"), emit: normal_germline_gvcf_index

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /usr/GenomeAnalysisTK.jar \
    -T HaplotypeCaller  \
    -R ${params.ref_fa} \
    -I ${normal_bam} \
    -o ${sampleID}_variants_raw.gvcf \
    -L ${params.target_gatk} \
    -stand_call_conf ${params.call_val} \
    -ERC GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000
    """
}