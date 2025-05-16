process QUALITY_CHECKS {
    tag "$sampleID"

    cpus 2
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.fragment_length_count.txt", mode: 'copy'
    
    input:
    tuple val(sampleID), file(sort_rm_filter_bam)

    output:
    tuple val(sampleID), file("*.fragment_length_count.txt")

    script:
    // Get the fragment length count from bam file for Quality Checks.
    """
    samtools view \
    -@ $task.cpus ${sort_rm_filter_bam[0]} \
    | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
    | sed -e 's/^[ \\t]*//' | awk -v sample="${sampleID}" -F' ' '{print sample,\$1,\$2}' OFS="\\t" > ${sampleID}.fragment_length_count.txt
    """
}
