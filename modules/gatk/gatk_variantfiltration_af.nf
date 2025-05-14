
process GATK_VARIANTFILTRATION_AF {
    // This modules is a port of the NYGC germline filtering scheme found at this site:
    // https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/germline/germline.wdl?at=7.4.0

    tag "$sampleID"

    cpus = 1
    memory = 10.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(idx)

    output:
    tuple val(sampleID), file("*haplotypecaller.gatk.af-gq-filtered.vcf"), emit: vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    ## Annotate FORMAT/AF
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp"  VariantAnnotator \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -O ${sampleID}_haplotypecaller.gatk.af.vcf.gz \
    -A AlleleFraction
    
    ## remove biallellic sites
    zcat ${sampleID}_haplotypecaller.gatk.af.vcf.gz \
    | awk '(\$5 !~ ",")' \
    > ${sampleID}.biallellic.vcf

    ## Variant filtration
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" VariantFiltration \
    -R ${params.ref_fa} \
    -V ${sampleID}.biallellic.vcf \
    -O ${sampleID}.haplotypecaller.af-gq-filtered.vcf.gz \
    --genotype-filter-name "AlleleFraction" \
    --genotype-filter-expression "(AF < 0.25 && AF > 0.0) || AF > 0.75" \
    --genotype-filter-name "GQ20" \
    --genotype-filter-expression "GQ < 20"

    ## filter with AF (deliver)
    zcat ${sampleID}.haplotypecaller.af-gq-filtered.vcf.gz \
    | grep -v "AlleleFraction" \
    > ${sampleID}_haplotypecaller.gatk.af-gq-filtered.vcf
    """
}
