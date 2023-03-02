process FUSIONCATCHER {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_4'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions': 'fusioncatcher' }", pattern: "*.{tsv,txt}", mode:'copy'


    input:
        tuple val(sampleID), path(reads)

    output:
        tuple val(sampleID), path("*_fusioncatcher_fusions.txt"), optional:true, emit: fusioncatcher_fusions
        tuple val(sampleID), path("*_fusioncatcher_summary.txt"), optional:true, emit: fusioncatcher_summary
        tuple val(sampleID), path("*_fusioncatcher.log"), emit: fusioncatcher_log

    script:

    def input_reads = reads.toString().replace(" ", ",")

    """
    fusioncatcher.py \\
        -d ${params.fusioncatcher_ref} \\
        -i ${input_reads} \\
        -p ${task.cpus} \\
        -o . \\
        --skip-blat \\
        --limitSjdbInsertNsj ${params.fusioncatcher_limitSjdbInsertNsj}

    mv final-list_candidate-fusion-genes.txt ${sampleID}_fusioncatcher_fusions.txt
    mv summary_candidate_fusions.txt ${sampleID}_fusioncatcher_summary.txt
    mv fusioncatcher.log ${sampleID}_fusioncatcher.log

    """



}
