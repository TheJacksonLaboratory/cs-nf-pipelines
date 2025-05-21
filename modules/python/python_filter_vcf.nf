process FILTER_VCF {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

    script:
    """
    python \
    ${projectDir}/bin/pta/vcf_filter.py \
    ${params.germline_filtering_vcf} \
    ${vcf} \
    ${sampleID}_final_filtered_${chrom}.vcf
    """
}
// NOTE: There are two similarly named scripts: vcf_filter.py and filter_vcf.py. 
//       The above script is used here. filter_vcf.py is used in gridss_chrom_filter.nf
