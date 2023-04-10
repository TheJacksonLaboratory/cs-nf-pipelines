process NANOSV {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"

    container 'quay.io/jaxcompsci/nanosv:1.2.4'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.nanosv_calls.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(bam), file(index)
    output:
        tuple val(sampleID), file("${sampleID}.nanosv_calls.vcf"), emit: nanosv_vcf
    script:
        """
        NanoSV -o ${sampleID}.nanosv_calls.vcf -t ${task.cpus} -c ${projectDir}/config/nanosv_config.ini ${bam}
        """
}