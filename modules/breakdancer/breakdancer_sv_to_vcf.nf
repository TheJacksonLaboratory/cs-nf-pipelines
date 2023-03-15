process BREAKDANCER_SV_TO_VCF {
    tag "$sampleID"
    
    cpus = 1
    memory = 20.GB
    time = "2:00:00"

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'
    
    input:
        tuple val(sampleID), file(bam), file(bai), file(breakdancer_sv)
    
    output:
        tuple val(sampleID), file("${sampleID}_breakdancer_sorted.vcf"), emit: breakdancer_vcf

    script:
        """
        python ${projectDir}/bin/breakdancer2vcfHeader.py -i ${breakdancer_sv} -o ${sampleID}_breakdancer.vcf
        bash ${projectDir}/bin/vcfSort.sh ${sampleID}_breakdancer.vcf ${sampleID}_breakdancer_sorted.vcf
        """
}