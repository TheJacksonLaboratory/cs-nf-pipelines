process REHEADER_VCR {
    tag "$sampleID"

    cpus = 1
    memory = 1.GB
    container "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_3"

    publishDir "${params.outdir}/unmerged_calls", pattern: "${sampleID}_${caller}Sort.vcf", mode: 'copy'

    input:
        tuple val(sampleID), file(vcf)
        val(caller)

    output:
        file("${sampleID}_${caller}Sort.vcf"), emit: vcf_rehead

    script:
        """
        printf "${sampleID}_${caller}\n" > rehead.txt
        bcftools reheader --samples rehead.txt \
            -o ${sampleID}_${caller}Sort.vcf \
            ${vcf}
        """	
}