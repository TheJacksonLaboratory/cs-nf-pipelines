process MAKE_GENOME_FILTER {
    tag "$sampleID"
    
    cpus 1
    memory 15.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ubuntu:20.04'

    publishDir "${params.pubdir}/genome_info", mode: 'copy'

    input:
    tuple val(sampleID), path(fai)
    file(blacklist)

    output:
    path('*.bed'), emit: bed
    path('*.sizes'), emit: sizes

    script: 
    fasta="\$(echo ${fai} | sed 's/.fai//g')"
    blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    """
    cut -f 1,2 ${fai} > ${fasta}.sizes
    $blacklist_filter > ${fasta}.include_regions.bed
    """
}
