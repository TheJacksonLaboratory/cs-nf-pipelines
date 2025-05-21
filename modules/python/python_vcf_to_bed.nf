process VCF_TO_BED {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.bed"), val(meta), val(chrom), emit: bed

    script:
    """
    python \
    ${projectDir}/bin/pta/vcf_to_bed.py \
    ${vcf} \
    | bedtools \
    merge \
    > ${sampleID}_candidate_merged_${chrom}.bed 
    """
}
