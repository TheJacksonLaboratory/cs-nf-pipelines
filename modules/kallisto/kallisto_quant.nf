process KALLISTO_QUANT {

    tag "$sampleID"

    cpus 12
    memory 84.GB
    time 24.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/kallisto:0.48.0--h15996b6_2'

    input:
        tuple val(sampleID), path(reads)

    output:
        tuple val(sampleID), path("*kallisto_quant.fusions.txt"), emit: kallisto_fusions
        tuple val(sampleID), path("*abundance.h5"), emit: kallisto_abundance

    script:
    """
    kallisto quant \
        -t $task.cpus \
        -i ${params.kallisto_index} \
        --fusion \
        -o . \
        ${reads}
    mv fusion.txt ${sampleID}.kallisto_quant.fusions.txt
    mv abundance.h5 ${sampleID}.abundance.h5
    """
}
// NOTE: 
// Index built with command: 
// singularity run /projects/omics_share/meta/containers/quay.io-biocontainers-kallisto-0.48.0--h15996b6_2.img kallisto index -k 31 -i Homo_sapiens.GRCh38.102.cdna.all.kallisto-0.48.0.index Homo_sapiens.GRCh38.102.cdna.all.fa.gz
