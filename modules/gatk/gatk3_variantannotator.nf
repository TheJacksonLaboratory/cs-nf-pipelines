process GATK_VARIANTANNOTATOR {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // Legacy Reasons Leave as GATK3 (public)
    // Flag --snpEffFile was removed in GATK4
    container 'broadinstitute/gatk3:3.6-0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.gen_org=='mouse' ? true : params.keep_intermediate

    input:
    tuple val(sampleID), file(sample_vcf), file(snpeff_vcf)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /usr/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R ${params.ref_fa} \
    -A SnpEff \
    --variant ${sample_vcf} \
    --snpEffFile ${snpeff_vcf} \
    -L ${sample_vcf} \
    -o ${sampleID}_GATKannotated.vcf
    """
}