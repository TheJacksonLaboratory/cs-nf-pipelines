process BCFTOOLS_GERMLINE_FILTER {
    // This modules is a port of the NYGC germline filtering scheme found at this site:
    // https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/germline/germline.wdl?at=7.4.0
    tag "$sampleID"
    
    cpus = 1
    memory = 2.GB
    time = '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*haplotypecaller.gatk.filtered.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*haplotypecaller.gatk.filtered.vcf.gz"), file("*haplotypecaller.gatk.filtered.vcf.gz.tbi"), emit: vcf_idx


    // NOTE: These are hard coded to resources provided at: https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/config/fasta_references.json
    //       Many of the files used in the filtering here are used again by VEP. Therefore, the reosource sets were combined to reduce the number of params. 

    script:
    """
    bgzip ${vcf}
    tabix -p vcf ${vcf}.gz

    ## Remove existing AF annotations from merged VCF
    bcftools annotate \
    -x INFO/AF \
    -Oz \
    ${vcf}.gz \
    > noaf.vcf.gz
    
    tabix -p vcf noaf.vcf.gz

    ## Annotate with NYGC AF for filtering
    bcftools annotate \
    --annotations ${params.vep_cache_directory}/annotations/04142020_NYGC_samples.vcf.gz \
    --columns 'INFO/AF,INFO/AC_Hom' \
    -Oz \
    noaf.vcf.gz \
    > ${sampleID}.final.annotated.vcf.gz

    tabix -p vcf ${sampleID}.final.annotated.vcf.gz

    ## filter variants >3% AF and >10 Homozygotes in NYGC vars
    bcftools filter \
    --exclude 'INFO/AF[*] > 0.03 || INFO/AC_Hom[*] > 10' \
    ${sampleID}.final.annotated.vcf.gz \
    > ${sampleID}.pop.filtered.vcf

    bgzip ${sampleID}.pop.filtered.vcf
    tabix -p vcf ${sampleID}.pop.filtered.vcf.gz

    ## select whitelist variants
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/vep_whitelist_38.20201118.vcf.gz \
    ${vcf}.gz \
    > ${sampleID}.whitelist.filtered.vcf.gz

    tabix -p vcf ${sampleID}.whitelist.filtered.vcf.gz

    ## select pgx variants
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/pgx_vep_hg38.vcf.gz \
    ${vcf}.gz \
    > ${sampleID}.pgx.filtered.vcf.gz

    tabix -p vcf ${sampleID}.pgx.filtered.vcf.gz

    ## select chd whitelist variants
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/chd_whitelist.vcf.gz \
    ${vcf}.gz \
    > ${sampleID}.chdwhitelist.filtered.vcf.gz

    tabix -p vcf ${sampleID}.chdwhitelist.filtered.vcf.gz

    ## select rwgs pgx variants
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/rWGS_PGx.bed.gz \
    ${vcf}.gz \
    > ${sampleID}.rwgspgx.filtered.vcf.gz

    tabix -p vcf ${sampleID}.rwgspgx.filtered.vcf.gz

    ## Select deep intronics
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/deep_intronic_whitelist_08132020.vcf.gz \
    ${vcf}.gz \
    > ${sampleID}.deep_intronics.filtered.vcf.gz

    tabix -p vcf ${sampleID}.deep_intronics.filtered.vcf.gz

    ## Select clinvar intronics
    bcftools view \
    -Oz \
    -R ${params.vep_cache_directory}/annotations/clinvar_deep_intronics_09012020.vcf.gz \
    ${vcf}.gz \
    > ${sampleID}.clinvar_intronics.filtered.vcf.gz

    tabix -p vcf ${sampleID}.clinvar_intronics.filtered.vcf.gz

    bcftools query -l ${vcf}.gz > samples.txt

    ## merge all filtered files for further processing
    bcftools concat \
    -a \
    -d all \
    ${sampleID}.pop.filtered.vcf.gz \
    ${sampleID}.whitelist.filtered.vcf.gz \
    ${sampleID}.pgx.filtered.vcf.gz \
    ${sampleID}.chdwhitelist.filtered.vcf.gz \
    ${sampleID}.rwgspgx.filtered.vcf.gz \
    ${sampleID}.deep_intronics.filtered.vcf.gz \
    ${sampleID}.clinvar_intronics.filtered.vcf.gz \
    | \
    bcftools view \
    -i 'GT[@samples.txt]="alt"' \
    | \
    bcftools sort \
    -Oz \
    > ${sampleID}_haplotypecaller.gatk.filtered.vcf.gz

    tabix -p vcf ${sampleID}_haplotypecaller.gatk.filtered.vcf.gz

    """
}