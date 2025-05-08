process GATK_FILTERMUECTCALLS {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*_mutect2_somatic.filtered.vcf.gz", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.filteringStats.tsv", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi), path(stats), path(contam_table), path(segments), file(read_orientation_model) //note: file() is used here on purpose. This is an optional input. 

    output:
    tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz"), file("*_mutect2_somatic.filtered.vcf.gz.tbi"), emit: mutect2_vcf_tbi
    tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    read_model_setting = params.ffpe ? "--orientation-bias-artifact-priors ${read_orientation_model}" : ""

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" FilterMutectCalls \
    -R ${params.ref_fa} \
    -V ${vcf} \
    --stats ${stats} \
    --contamination-table ${contam_table} \
    --tumor-segmentation ${segments} \
    ${read_model_setting} \
    -O ${sampleID}_mutect2_somatic.filtered.vcf.gz
    """
}
