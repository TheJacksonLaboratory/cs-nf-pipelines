process SNPEFF_ONEPERLINE {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '00:10:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    input:
    tuple val(sampleID), file(vcf)
    val(indel_snp)

    output:
    tuple	val(sampleID), file("*.vcf"), emit: vcf

    script:
    if (indel_snp == 'INDEL'){
        output_suffix = 'INDEL_snpeff.vcf'
    }
    if (indel_snp =='SNP'){
        output_suffix = 'SNP_snpeff.vcf'
    }
    if (indel_snp == 'BOTH'){
        output_suffix = 'snp_indel_snpeff.vcf'
    }
    """
    cat ${vcf} | perl ${projectDir}/bin/shared/vcfEffOnePerLine.pl > ${sampleID}_oneperline_${output_suffix}
    """
}
