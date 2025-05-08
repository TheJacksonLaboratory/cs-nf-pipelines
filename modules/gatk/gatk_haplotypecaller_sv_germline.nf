process GATK_HAPLOTYPECALLER_SV_GERMLINE {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.*vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(read_name), path(interval), val(index)
    
    output:
    tuple val(sampleID), path("*.*vcf"), emit: vcf
    tuple val(sampleID), path("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" HaplotypeCaller  \
    -R ${params.ref_fa} \
    -I ${normal_bam} \
    -O ${sampleID}_${index}_variants_raw.gvcf \
    -L ${interval} \
    -XL ${params.excludeIntervalList} \
    -stand-call-conf ${params.call_val} \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    -G AS_StandardAnnotation \
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
    -ERC GVCF
    """
}
