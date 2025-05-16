process GATKv3_5_GENOTYPEGVCF {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk3:3.5-0'
    
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(normal_germline_gvcf)
    tuple val(sampleID), file(normal_germline_gvcf_index)

    output:
    tuple val(sampleID), file("*.*vcf"), emit: normal_germline_vcf
    tuple val(sampleID), file("*.vcf.idx"), emit: normal_germline_vcf_index

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    
    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R ${params.ref_fa} \
    --variant ${sampleID}_variants_raw.gvcf \
    -o ${sampleID}_variants_raw.vcf
    """
    }
   