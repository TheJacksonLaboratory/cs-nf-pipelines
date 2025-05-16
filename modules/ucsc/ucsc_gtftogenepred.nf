process UCSC_GTFTOGENEPRED {
    cpus 1
    memory 24.GB
    time 1.h
    
    tag "$gtf"

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/ucsc-gtftogenepred:447--h954228d_0"

    publishDir "${params.pubdir}", pattern: "${gtf.baseName}.refFlat.txt", mode:'copy'

    input:
        path(gtf)

    output:
        path("${gtf.baseName}.refFlat.txt"), emit: refFlat

    script:
        """
        gtfToGenePred \
            -genePredExt \
            ${gtf} \
            -ignoreGroupsWithoutExons \
            tmp_output.txt

        cat tmp_output.txt | awk -v OFS="\t" '{print \$12,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' > ${gtf.baseName}.refFlat.txt
        """
}

// Module from: https://github.com/KU-GDSC/workflows
