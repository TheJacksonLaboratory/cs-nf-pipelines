process GATKv3_5_ApplyRecalibration {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk3:3.5-0'
    
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(normal_germline_vcf)
    tuple val(sampleID), file(normal_germline_vcf_index)
    tuple val(sampleID), file(normal_germline_recal)
    tuple val(sampleID), file(normal_germline_tranches)

    output:
    tuple val(sampleID), file("*.*recalibrated.filtered.vcf"), emit: normal_germline_recalibrated_vcf
    tuple val(sampleID), file("*.*recalibrated.filtered.vcf.idx"), emit: normal_germline_recalibrated_vcf_index

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    
    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R ${params.ref_fa} \
    -input ${sampleID}_variants_raw.vcf \
    --ts_filter_level 99.6 \
    -tranchesFile ${sampleID}.tranches.txt \
    -recalFile ${sampleID}.recal.txt \
    -mode SNP
    -o ${sampleID}_variants_raw.recalibrated.filtered.vcf
    """
}