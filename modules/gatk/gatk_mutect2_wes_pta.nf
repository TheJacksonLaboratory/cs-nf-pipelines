process GATK_MUTECT2 {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), file("*_somatic.vcf.gz"), file("*_somatic.vcf.gz.tbi"), file("*.stats"), emit: vcf_tbi_stats
    tuple val(sampleID), file("*f1r2.tar.gz"), emit: f1r2

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    germline_genotype = params.genotype_germline ? '--genotype-germline-sites true' : ''
    pon_genotype = params.genotype_pon ? '--genotype-pon-sites true' : ''

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" Mutect2 \
    -R ${params.ref_fa} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -normal ${normal_name} \
    --germline-resource ${params.gnomad_ref} \
    --panel-of-normals ${params.pon_ref} \
    --f1r2-tar-gz ${sampleID}.f1r2.tar.gz \
    ${germline_genotype} \
    ${pon_genotype} \
    --pileup-detection \
    --dont-use-soft-clipped-bases false \
    -L ${params.target_gatk} \
    --native-pair-hmm-threads 4 \
    --annotation QualByDepth \
    --annotation RMSMappingQuality \
    --annotation FisherStrand \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest \
    --min-base-quality-score 20 \
    -O ${sampleID}_mutect2_somatic.vcf.gz
    """
}

/*
As of v4.1, there is no longer a need to specify the tumor sample name with -tumor. You need only specify the normal sample name with -normal, if you include a normal.

Starting with v4.0.4.0, GATK recommends the default setting of --af-of-alleles-not-in-resource, which the tool dynamically adjusts for different modes. 
tumor-only calling sets the default to 5e-8, tumor-normal calling sets it to 1e-6 and mitochondrial mode sets it to 4e-3. 
For previous versions, the default was 0.001, the average heterozygosity of humans. 
For other organisms, change --af-of-alleles-not-in-resource to 1/(ploidy*samples in resource).

https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

*/
