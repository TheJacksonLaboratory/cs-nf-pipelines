process DELLY_POST_PROCESS {
    tag "$sampleID"
    
    cpus = 1
    memory = 20.GB
    time = "2:00:00"

    container 'quay.io/biocontainers/bcftools:1.10.2--h4f4756c_3'
    
    input:
        tuple val(sampleID), file(delly_bcf)
    
    output:
        tuple val(sampleID), file("${sampleID}_delly.vcf"), emit: delly_vcf

    script:
        """
        bcftools view ${delly_bcf} > "${sampleID}_delly_unsorted.vcf"
        bash ${projectDir}/bin/vcfSort.sh "${sampleID}_delly_unsorted.vcf" ${sampleID}_delly.vcf
        """
}