process GATK_FILTER_VARIANT_TRANCHES {
    // This modules is a port of the NYGC germline filtering scheme found at this site:
    // https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/germline/germline.wdl?at=7.4.0
    
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '01:30:00'
    errorStrategy 'ignore'

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*vcf", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)
    tuple val(sampleID), file(vcf_index)

    output:
    tuple val(sampleID), file("*.*vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    gatk --java-options "-Xmx${my_mem}G" FilterVariantTranches  \
    -V ${vcf} \
    -O ${sampleID}_haplotypecaller.gatk.filtered.genotypedGVCFs.vcf \
    -L ${params.target_gatk} \
    --snp-tranche 99.9 --snp-tranche 99.95 \
    --indel-tranche 99.0 --indel-tranche 99.4 \
    --resource ${params.hapmap} \
    --resource ${params.omni} \
    --resource ${params.onekG} \
    --resource ${params.dbsnp} \
    --info-key CNN_1 \
    --create-output-variant-index true
    """
}