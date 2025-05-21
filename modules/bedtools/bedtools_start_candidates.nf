process BEDTOOLS_STARTCANDIDATES {
    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

    script:

    """
    bedtools \
    intersect \
    -header \
    -a ${vcf} \
    -b ${params.intervalListBed} \
    -v \
    > ${sampleID}_startCand_merged_${chrom}.vcf
    """
}
