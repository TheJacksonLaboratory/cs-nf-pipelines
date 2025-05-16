process GERMLINE_VCF_FINALIZATION {
    tag "$sampleID"

    cpus 1
    memory 5.GB
    time 1.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*final.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)
    val(filtered)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf

    script:

    output_suffix = filtered == 'filtered' ?  'filtered' : 'unfiltered'

    """
    python \
    ${projectDir}/bin/pta/annotate_id.py \
    ${vcf} \
    ${sampleID}_germline_vep_cosmic_cancerResitMut_annotated_id.vcf

    python \
    ${projectDir}/bin/pta/rename_csq_vcf.py \
    ${sampleID}_germline_vep_cosmic_cancerResitMut_annotated_id.vcf \
    ${sampleID}_germline_snv_indel_annotated_${output_suffix}_final.vcf
    """
}
