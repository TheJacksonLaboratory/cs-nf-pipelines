process GENERATE_RRNA_INTERVALS {
    tag "$gff"

    cpus 1
    memory 5.GB 
    time '00:20:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/perl:0.1.0"

    publishDir "${params.pubdir}", pattern: "${gff.baseName}.rRNA_intervals.list", mode:'copy'

    input:
        path(gff)
        path(dict)
    output:
        path("${gff.baseName}.rRNA_intervals.list"), emit: rRNA_intervals

    script:
        """
        grep "@" ${dict} > ${gff.baseName}.rRNA_intervals.list

        grep 'gene_biotype' ${gff} | \
            grep "rRNA" | \
            sed -e 's/;/\t/g' | \
            cut -f 1,4,5,7,9 | \
            sed 's/ID=//g' >> ${gff.baseName}.rRNA_intervals.list
        """
}

// Module from: https://github.com/KU-GDSC/workflows
