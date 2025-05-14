process GRIPSS_SOMATIC_FILTER {
    tag "$sampleID"

    cpus = 1
    memory = 30.GB
    time = '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/hmftools-gripss:2.3.2--hdfd78af_0'

    stageInMode = 'copy'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*gripss.filtered.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), val(meta), val(normal_name), val(tumor_name)


    output:
    tuple val(sampleID), path('*gripss.filtered.vcf.gz'),path('*gripss.filtered.vcf.gz.tbi'), val(meta), val(normal_name), val(tumor_name), val('gridss'), emit: gripss_filtered_bgz 
                                                                                                                                            //note: while this is "GRIPSS" filtering.  
                                                                                                                                            //      GRIDSS was the caller and downstream 
                                                                                                                                            //      scripts expect "gridss" as the tool name. 
    tuple val(sampleID), path('*.gripss.vcf.gz'), path('*.gripss.vcf.gz.tbi'), val(meta), val(normal_name), val(tumor_name), val('gridss'), emit: gripss_all_bgz

    script:
    """
    gripss -Xmx29g \
        -sample ${tumor_name} \
        -reference ${normal_name} \
        -ref_genome_version 38 \
        -ref_genome ${params.ref_fa} \
        -pon_sgl_file ${params.gripss_pon}/sgl_pon.38.bed \
        -pon_sv_file ${params.gripss_pon}/sv_pon.38.bedpe \
        -known_hotspot_file ${params.gripss_pon}/known_fusions.38.bedpe \
        -repeat_mask_file ${params.gripss_pon}/repeat_mask_data.38.fa.gz \
        -vcf ${vcf} \
        -output_dir .

        mv ${tumor_name}.gripss.filtered.vcf.gz ${sampleID}.gripss.filtered.vcf.gz
        mv ${tumor_name}.gripss.filtered.vcf.gz.tbi ${sampleID}.gripss.filtered.vcf.gz.tbi
        mv ${tumor_name}.gripss.vcf.gz ${sampleID}.gripss.vcf.gz
        mv ${tumor_name}.gripss.vcf.gz.tbi ${sampleID}.gripss.vcf.gz.tbi

    """
}
