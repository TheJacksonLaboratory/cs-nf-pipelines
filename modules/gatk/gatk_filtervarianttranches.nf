process GATK_FILTER_VARIANT_TRANCHES {
    // This modules is a port of the NYGC germline filtering scheme found at this site:
    // https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/germline/germline.wdl?at=7.4.0
    
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '01:30:00'
    errorStrategy 'ignore'

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.*vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(vcf_index)

    output:
    tuple val(sampleID), file("*.*vcf"), file("*.idx"), emit: vcf_idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" FilterVariantTranches  \
    -V ${vcf} \
    -O ${sampleID}_haplotypecaller.gatk.filtered.genotypedGVCFs.vcf \
    --snp-tranche 99.9 --snp-tranche 99.95 \
    --indel-tranche 99.0 --indel-tranche 99.4 \
    --resource ${params.hapmap} \
    --resource ${params.omni} \
    --resource ${params.phase1_1000G} \
    --resource ${params.dbSNP} \
    --info-key CNN_1D \
    --create-output-variant-index true
    """
}
