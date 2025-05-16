process FUSION_REPORT {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time 2.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/fusion-report:2.1.5--pyhdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/fusion-report/'}", mode:'copy'

    input:
        tuple val(sampleID), path(arriba), path(fusioncatcher), path(jaffa), path(pizzly), path(squid), path(starfusion)

    output:
        tuple val(sampleID), file("${sampleID}_fusion_list.tsv"), emit: fusion_inspector_input_list
        tuple val(sampleID), file("${sampleID}_fusion_genes_mqc.json"), emit: summary_fusions_mq
        tuple val(sampleID), file("*"), emit: report
    
    script:
    def extra_params = params.fusion_report_opt ? params.fusion_report_opt : ''
    def tools =  !arriba.empty() ? "--arriba ${arriba} " : ''
        tools += !jaffa.empty() ? "--jaffa ${jaffa} " : ''
        tools += !fusioncatcher.empty() ? "--fusioncatcher ${fusioncatcher} " : ''
        tools += !pizzly.empty() ? "--pizzly ${pizzly} " : ''
        tools += !squid.empty() ? "--squid ${squid} " : ''
        tools += !starfusion.empty() ? "--starfusion ${starfusion} " : ''

    """
    fusion_report run ${sampleID} . ${params.databases} ${tools} ${extra_params}
    mv fusion_list.tsv ${sampleID}_fusion_list.tsv
    mv fusion_list_filtered.tsv ${sampleID}_fusion_list_filtered.tsv
    mv fusion_genes_mqc.json ${sampleID}_fusion_genes_mqc.json
    """
}
