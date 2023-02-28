process KALLISTO_QUANT {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/biocontainers/kallisto:0.48.0--h15996b6_2'

    input:
        tuple val(sampleID), file(reads)

    output:
        tuple val(sampleID), file("${sampleID}.kallisto_quant.fusions.txt"), emit: kallisto_fusions

    script:
    """
    kallisto quant \
        -t $task.cpus \
        -i ${params.kallisto_index} \
        --fusion \
        -o . \
        ${reads}
    mv fusion.txt ${sampleID}.kallisto_quant.fusions.txt
    """
}
// NOTE: 
// Index built with command: 
// singularity run /projects/omics_share/meta/containers/quay.io-biocontainers-kallisto-0.48.0--h15996b6_2.img kallisto index -k 31 -i Homo_sapiens.GRCh38.102.cdna.all.kallisto-0.48.0.index Homo_sapiens.GRCh38.102.cdna.all.fa.gz
