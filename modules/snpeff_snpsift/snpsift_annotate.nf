process SNPSIFT_ANNOTATE {
    tag "$sampleID"

    cpus = 1
    memory = 20.GB
    time = '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'
    
    publishDir "${params.pubdir}/${sampleID}", pattern:"*dbsnpID.vcf", mode:'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern:"*_germline_snv_indel_annotated_filtered_final.vcf", mode:'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern:"*.vcf", mode:'copy', enabled: params.workflow == 'amplicon' ? true : false

    input:
    tuple val(sampleID), file(vcf)
    path(annot_source)
    path(annot_index)
    val(output_suffix)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:

    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    output_name = params.gen_org == 'mouse' && params.workflow =='pta' ? "${sampleID}_germline_snv_indel_annotated_filtered_final.vcf" : "${vcf.baseName}_${output_suffix}.vcf"

    """
    java -Xmx${my_mem}G -jar /opt/snpEff/SnpSift.jar \
    annotate -noDownload -id ${annot_source} ${vcf} > ${output_name}
    """
}
