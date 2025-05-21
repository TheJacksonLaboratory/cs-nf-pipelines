process MINIMAP2_MAP_ONT {
    tag "$sampleID"
    
    cpus 16
    memory 24.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io-biocontainers-minimap2-2.24--h7132678_1'

    publishDir "${params.pubdir}/${sampleID + '/alignments'}", pattern: "*_aln.sam", mode:'copy', enabled: params.keep_intermediate ? true : false

    input:
        tuple val(sampleID), file(fastq)
        path(minimap2_index)
    output:
        tuple val(sampleID), file("${sampleID}_aln.sam"), emit: minimap_sam
    script:
        """
        minimap2 --secondary=no -t ${task.cpus} --MD -ax map-ont ${minimap2_index} ${fastq}  > ${sampleID}_aln.sam
        """
}
