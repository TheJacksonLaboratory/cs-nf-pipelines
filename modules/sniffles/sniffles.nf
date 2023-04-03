process SNIFFLES {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"

    container 'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0'

    publishDir "${params.pubdir}/unmerged_calls", pattern: "${sampleID}.sniffles_calls.vcf", mode: "copy"

    input:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai)
    output:
        tuple val(sampleID), file("${sampleID}.sniffles_calls.vcf"), emit: sniffles_vcf
    script:
        """
        sniffles --input ${pbmm2_bam} --vcf ${sampleID}.sniffles_calls.vcf
        """
}